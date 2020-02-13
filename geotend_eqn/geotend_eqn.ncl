load "./src_test/poisson_solver.ncl"

; Phi, J, T

y4m2 = "201607"
path = "/data2/eraint/1.5x1.5/"

; Read file
; fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
; U  = fU->U_GDS0_ISBL 
; fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
; V  = fV->V_GDS0_ISBL 
fin  = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T    = fin->T_GDS0_ISBL(1,{200:1000},:,:)    ; Temperature [K]
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Z    = fin->Z_GDS0_ISBL(1,{200:1000},:,:)    ; Geopotential
fin  = addfile(path+"PHYpl/"+"Tphy.pl."+y4m2+".nc", "r")
J_cp = fin->Tphy(1,{200:1000},:,:)/(3600.*6) ; =J/Cp [K/6hour] -> [K/s]

lon = T&g0_lon_3
lat = T&g0_lat_2 
prs = T&lv_ISBL1    ; [hPa] -> [Pa]
a   = 6.37122e06    ; [m]
R   = 287.          ; 
pi  = 4*atan(1.)
d2r = pi/180.       ; Degree to radian
OMG = 7.2921e-5
P   = conform(T, prs, 0)

; Calculate Static stability 
; make only zonal average for practical reason
Tzm = dim_avg(T)
sig = static_stability(prs, Tzm, 0, 0)
sig = (R/conform(sig, prs, 0))*sig  ; sig(t,p,y)

; Calculate f
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat
f({0}) = 1.e-5          ; practical reason
                        ; Assume! f_equator is small value.
f3   = conform(T, f, 1)
beta = 2*OMG*cos(lat*d2r)/a

; Calculate geostrophic wind Ug, Vg
; Use f instead f0.
Ug = -(1./f3)*center_finite_diff_n(Z, a*lat*d2r, False, 1, 1)
Vg =  (1./f3)*center_finite_diff_n(Z, a*lon*d2r, False, 1, 2)
Vg = Vg/conform(Vg, cos(lat*d2r), 1)

; Calculate geostrophic vorticity
; vo = uv2vr_cfd(Ug, Vg, lat, lon, 3)
vo = uv2vrF(Ug, Vg)

; Forcing 1 : Absolute vorticity advection
grid = 1       ; 1 : fixed grid, 0 : gaussian grid
avo  = f3+vo   ; Absolute vorcity
copy_VarCoords(T, avo)
gridType = 1   ; 1 : fixed grid
longName = "Absolute vorticity advection"
units = "1/s^2"
if (lat(0) .gt. lat(1)) then
    avo_adv  = advect_variable(Ug(:,::-1,:), Vg(:,::-1,:), \
                              avo(:,::-1,:), 1, longName, units, 0)
    avo_adv  = avo_adv(:,::-1,:)
else
    avo_adv = advect_variable(Ug, Vg, avo, 1, longName, units, 0)
end if    
F1 = -f3*avo_adv

; Forcing 2 : Thickness advection
dZdp = R*T/P  ; -dgeo/dp = R*T/P
copy_VarCoords(T, dZdp)
sig  = where(sig .le. 0, 1.e-7, sig)   ; Assume! remove negative sigma.
sig3 = conform_dims(dimsizes(T), sig, (/0,1/))
if (lat(0) .gt. lat(1)) then
    dZdp_adv = advect_variable(Ug(:,::-1,:), Vg(:,::-1,:), \
                             dZdp(:,::-1,:), 1, longName, units, 0)
    dZdp_adv = dZdp_adv(:,::-1,:)
else
    dZdp_adv = advect_variable(Ug, Vg, dZdp, 1, longName, units, 0)
end if    
F2   = (-(f3^2.)/sig3)*dZdp_adv
; optional (diff p)
F2   = -center_finite_diff_n(F2, prs, False, 1, 0)
; log-P
; logp = log(prs)       ; Finite difference in log-P (d/dp -> 1/p*d/dlogp.)
; F2   = center_finite_diff_n(F2, logp, False, 1, 1)
; F2   = -F2/P          ; take it back to pressure

; Forcing 3 : Heating
F3   = ((R/P)*f3^2./sig3*J_cp) 
; optional (diff p)
F3   = -center_finite_diff_n(F3, prs, False, 1, 0)
; log-P
; F3   = center_finite_diff_n(F3, logp, False, 1, 1)
; F3   = -F3/conform(F3, prs, 1)          ; take it back to pressure

; =============== Need to work================
;================================================================

; Tendency term X = (dZ/dt)
; Computes using lapsF function in ncl.
; For test:
;     ts    = shagC(tg)
;     n     = ispan(0,nlat-1,1)
;     nn1   = conform(ts,n*(n+1),1)
;     lap2s = -ts*nn1/(R^2)
;     lap2  = shsgC(lap2s,nlon)
; delZ = lapsF(Z)
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Z   := fin->Z_GDS0_ISBL(0:2,{200:1000},:,:)    ; Geopotential
X    = center_finite_diff_n(Z, 21600., False, 1, 0)
F = F1 + F2 + F3
copy_VarCoords(Z(1,:,:,:), F)

opt = True
opt@T0 = dim_avg_Wrap(dim_avg_Wrap(X(1,:,:,:)))
opt@lev_units = "hPa"
opt@lev_units = "hPa"
opt@dt = 2.e-1
opt@thres = 1.047e16
solution = poisson_solver(F, opt)
printVarSummary(solution)

; Solve equation using SOR(Successive Over Relexation) method

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
fc->Z   = Z
fc->solution = solution
fc->forcing = F
















