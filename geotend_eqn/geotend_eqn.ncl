load "./src_test/poisson_solver.ncl"

y4m2 = "201602"
path = "/data2/eraint/1.5x1.5/"
time_idx = 7

; Read file
; fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
; U  = fU->U_GDS0_ISBL 
; fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
; V  = fV->V_GDS0_ISBL 
fin  = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T    = fin->T_GDS0_ISBL(time_idx,{200:1000},:,:)    ; Temperature [K]
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
PHI  = fin->Z_GDS0_ISBL(time_idx,{200:1000},:,:)    ; Geopotential [m**2 s**-2]
fin  = addfile(path+"PHYpl/"+"Tphy.pl."+y4m2+".nc", "r")
J_cp = fin->Tphy(time_idx,{200:1000},:,:)/(3600.*6) ; =J/Cp [K/6hour] -> [K/s]
; test Q1, Q2
; cp   = 1005.       ; [J/kg/K]
; fin  = addfile("./apparent_heat_moist.nc", "r")
; J_cp = fin->Q1(time_idx,{200:1000},:,:)/cp  ; [J/kg/s] -> [K/s]

; load coordinate
lon = T&g0_lon_3
lat = T&g0_lat_2 
prs = (T&lv_ISBL1)*100.   ; [hPa] -> [Pa]
a   = 6.37122e06          ; [m]
R   = 287.          ; 
pi  = 4*atan(1.)
d2r = pi/180.       ; Degree to radian
OMG = 7.2921e-5
P   = conform(T, prs, 0)

; Calculate Static stability 
; make only zonal average for practical reason
Tzm = dim_avg(dim_avg(T))
sp0 = static_stability(prs, Tzm, 0, 0)
sig = (R/conform(sp0, prs, 0))*sp0  ; sig(t,p,y)

; Calculate f
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat                 ; practical reason
f({0}) = 1.e-5              ; Assume! f_equator is small value.
f3   = conform(T, f, 1)

; Forcing 1 : Absolute vorticity advection
gtype = 1   ; 1 = fixed grid (default), 0 = gaussian grid
vdiff = 1   ; 1 = pressure (default),   0 = logp,

; Calculate geostrophic wind Ug, Vg
; Use f instead f0.
Ug = -(1./f3)*center_finite_diff_n(PHI, a*lat*d2r, False, 1, 1)
Vg =  (1./f3)*center_finite_diff_n(PHI, a*lon*d2r, False, 1, 2)
Vg = Vg/conform(Vg, cos(lat*d2r), 1)

; Calculate geostrophic vorticity
if (gtype .eq. 1) then vo = uv2vrF(Ug, Vg)
else                   vo = uv2vrG(Ug, Vg) end if    

avo  = f3+vo   ; Absolute vorcity
copy_VarCoords(T, avo)
longName = "Absolute vorticity advection"
units = "1/s^2"
if (lat(0) .gt. lat(1)) then latord = -1
else                         latord =  1 end if
avo_adv = -advect_variable(Ug(:,::latord,:), \ ; advect_variable 
                           Vg(:,::latord,:), \ ; = u*d/dx + v*d/dy
                          avo(:,::latord,:), \ ; Need to add minus sign
                 gtype, longName, units, 0)
avo_adv = avo_adv(:,::latord,:)
F1      = f3*avo_adv

; Forcing 2 : Thickness advection 
dPHI_dp = -R*T/P    ; replace dphi_dp using temperature
copy_VarCoords(T, dPHI_dp)
sig  = where(sig .le. 0, 1.e-7, sig)   ; Assume! remove negative sigma.
sig3 = conform_dims(dimsizes(T), sig, (/0/))
_dPHI_dp = -dPHI_dp     ; Avoid for warning
copy_VarCoords(T, _dPHI_dp)
temp_adv = -advect_variable(Ug(:,::latord,:), \ ; advect_variable 
                            Vg(:,::latord,:), \ ; = u*d/dx + v*d/dy
                      _dPHI_dp(:,::latord,:), \ ; Need to add minus sign
                  gtype, longName, units, 0)
temp_adv = temp_adv(:,::latord,:)
F2 = ((f3^2.)/sig3)*temp_adv

if (vdiff .eq. 1) then  ; Diff in pressure
    F2 = -center_finite_diff_n(F2, prs, False, 1, 0)
else  ; Diff in log-P (d/dp -> 1/p*d/dlogp.)
    logp = log(prs)       
    F2 = -center_finite_diff_n(F2, logp, False, 1, 0)
    F2 = F2/P        ; take it back to pressure
end if    


; Forcing 3 : Heating
J_cp@_FillValue = 0.
delete(J_cp@_FillValue)
F3 = (R/P)*((f3^2.)/sig3)*J_cp
if (vdiff .eq. 1) then  ; Diff in pressure
    F3 = -center_finite_diff_n(F3, prs, False, 1, 0)
else  ; Diff in log-P (d/dp -> 1/p*d/dlogp.)
    F3 = -center_finite_diff_n(F3, logp, False, 1, 0)
    F3 = F3/P          ; take it back to pressure
end if
;================================================================
; Inverse laplacian
F = F1 + F2 + F3
; F = F3
copy_VarCoords(Zg, F)

bt = time_idx-1
ft = time_idx+1

; Tendency term X = (dZ/dt)
; Computes using lapsF function in ncl.
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Zg  := fin->Z_GDS0_ISBL(bt:ft,{200:1000},:,:)    ; Geopotential
Z   := Zg
copy_VarCoords(Zg, Z)
X    = center_finite_diff_n(Z, 21600., False, 1, 0)

; Solve equation using SOR(Successive Over Relexation) method
opt = True
opt@T0 = dim_avg_Wrap(dim_avg_Wrap(T(:,:,:)))
opt@lev_units = "hPa"
opt@dt = 1.e8
opt@f0 = 1.e-4
opt@thres = 1.e-10
opt@upper_bound = X(1,0,:,:)
opt@lower_bound = X(1,22,:,:)
solution = poisson_solver(F, opt)

; ============== test ==============
copy_VarCoords(T, F1)
copy_VarCoords(T, F2)
copy_VarCoords(T, F3)
copy_VarCoords(T, X(1,:,:,:))
ofile = "./geotend.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->F1  = F1
fc->F2  = F2
fc->F3  = F3
fc->X   = X(1,:,:,:)
fc->Z   = Z(1,:,:,:)
fc->solution = solution
fc->forcing = F
















