load "./src_test/poisson_solver.ncl"

; Calculate forcing in tendency equation /*{{{*/
undef("geotend_forcing")
function geotend_forcing(PHI[*][*][*], \ ; Geopotential [m2 s-2]
                        J_cp[*][*][*], \ ; DT/Dt ... = J/Cp [K/s]
                        opt:logical)
;
;  - Input argument
;    PHI  : Geopotential [lev] x [lat] x [lon]
;    J_cp : = J/Cp       [lev] x [lat] x [lon]
;    opt  : Option (See below)
;
;  - Description : Compute forcing in tendency equation
;    F1 : Absolute vorticity advection
;    F2 : Thickness advection 
;    F3 : Heating
;
;  - Option
;    opt@T         = T([lev]x[lat]x[lon], [K]) 
;       @sig_lat   = True: σ(z,y), False: σ(z) (default: False)
;       @lev_units = "hPa", "Pa"               (default: "Pa")
;       @f0        = Coriolis parameter
;       @gtype     = Grid type  
;                    1: fixed grid (default), 0: gaussian grid
;       @vdiff     = vertical differential
;                    1: pressure (default),   0: logp
;
begin

    ;-- Set variables
a   = 6.37122e06    ; Earth radius [m]
R   = 287.          ; Gas constant [J/kg/K] 
pi  = 4*atan(1.)    ; = 3.141592...
d2r = pi/180.       ; Degree to radian
OMG = 7.2921e-5

lev_units = "Pa"
sig_lat = False
gtype = 1   ; 1 = fixed grid (default), 0 = gaussian grid
vdiff = 1   ; 1 = pressure (default),   0 = logp,

if (isatt(opt,    "sig_lat")) then sig_lat   = opt@sig_lat end if
if (isatt(opt,      "gtype")) then gtype     = opt@gtype   end if
if (isatt(opt,      "vdiff")) then vdiff     = opt@vdiff   end if
if (isatt(opt,  "lev_units")) then lev_units = opt@lev_units 
else
    print("Note! opt@lev_units = Pa is assumed, change opt if not.")
end if

; load coordinate
dnam = getvardims(PHI)
zdim = dnam(0)
ydim = dnam(1)
xdim = dnam(2)
lon  = PHI&$xdim$
lat  = PHI&$ydim$
prs  = PHI&$zdim$
if (lev_units .eq. "hPa") then prs = prs*100. end if    
    
; Calculate Static stability 
logp = log(prs)       
T    = -(1/R)*center_finite_diff_n(PHI, logp, False, 1, 0)
if (isatt(opt, "T")) then T = opt@T end if

Tzm  = dim_avg(dim_avg(T))   ; make only zonal average
sig_dims = (/0/)             ; for practical reason
if (sig_lat) then Tzm      := dim_avg(T) end if
if (sig_lat) then sig_dims := (/0,1/)    end if

sp0  = static_stability(prs, Tzm, 0, 0)
sig  = (R/conform(sp0, prs, 0))*sp0    ; sig(t,p,y)
sig  = where(sig .le. 0, 1.e-7, sig)   ; Assume! remove negative sigma.
sig3 = conform_dims(dimsizes(PHI), sig, sig_dims)

; Calculate f (Default: Use f instead f0.)
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat                 ; practical reason
f({0}) = 1.e-5              ; Assume! f_equator is small value.
f3  = conform(PHI, f, 1)
if (isatt(opt, "f0")) then f3 = opt@f0 end if

; Calculate geostrophic wind Ug, Vg
Ug = -(1./f3)*center_finite_diff_n(PHI, a*lat*d2r, False, 1, 1)
Vg =  (1./f3)*center_finite_diff_n(PHI, a*lon*d2r, False, 1, 2)
Vg = Vg/conform(Vg, cos(lat*d2r), 1)

;-- Main: Calculate Forcing
; Forcing 1 : Absolute vorticity advection
; Calculate geostrophic vorticity
if (gtype .eq. 1) then vo = uv2vrF(Ug, Vg)
else                   vo = uv2vrG(Ug, Vg) end if    
avo  = f3+vo   ; Absolute vorcity
copy_VarCoords(PHI, avo)
longName = "Absolute vorticity advection"
units = "1/s^2"
if (lat(0) .gt. lat(1)) then latord = -1
else                         latord =  1 end if
avo_adv = -advect_variable(Ug(:,::latord,:), \ ; advect_variable 
                           Vg(:,::latord,:), \ ; = u*d/dx + v*d/dy
                          avo(:,::latord,:), \ ; Need to add minus sign
                 gtype, longName, units, 0)
avo_adv = avo_adv(:,::latord,:)
F1 = f3*avo_adv


; Forcing 2 : Thickness advection 
P        = conform(PHI, prs, 0)
dPHI_dp  = -R*T/P       ; replace dphi_dp using temperature
_dPHI_dp = -dPHI_dp     ; Avoid for warning
copy_VarCoords(PHI, _dPHI_dp)
temp_adv = -advect_variable(Ug(:,::latord,:), \ ; advect_variable 
                            Vg(:,::latord,:), \ ; = u*d/dx + v*d/dy
                      _dPHI_dp(:,::latord,:), \ ; Need to add minus sign
                  gtype, longName, units, 0)
temp_adv = temp_adv(:,::latord,:)
F2 = ((f3^2.)/sig3)*temp_adv

if (vdiff .eq. 1) then  ; Diff in pressure
    F2 = -center_finite_diff_n(F2,  prs, False, 1, 0)
else  ; Diff in log-P (d/dp -> 1/p*d/dlogp.)
    F2 = -center_finite_diff_n(F2, logp, False, 1, 0)
    F2 = F2/P        ; take it back to pressure
end if    


; Forcing 3 : Heating
J_cp@_FillValue = 0.
delete(J_cp@_FillValue)
F3 = (R/P)*((f3^2.)/sig3)*J_cp
if (vdiff .eq. 1) then  ; Diff in pressure
    F3 = -center_finite_diff_n(F3,  prs, False, 1, 0)
else  ; Diff in log-P (d/dp -> 1/p*d/dlogp.)
    F3 = -center_finite_diff_n(F3, logp, False, 1, 0)
    F3 = F3/P          ; take it back to pressure
end if

    F = F1 + F2 + F3
    copy_VarCoords(PHI, F)
    F@F1 = F1   ; Absolute vorticity advection
    F@F2 = F2   ; Thickness advection 
    F@F3 = F3   ; Heating

    return F
end ; /*}}}*/
