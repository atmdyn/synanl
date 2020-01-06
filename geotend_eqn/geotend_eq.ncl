; ; /*{{{*/
; undef("laplcaian")
; function laplcaian(Y)
; ; argv = [time] x [lev] x [lat] x [lon]
; begin
;     
;     dx  = Y&lon
;     Yxx = (Y(:,:,:,ix+1) - 2*Y(:,:,:,ix) + Y(:,:,:,ix-1)) / (dx**2)
;
;     return 0
; end ; /*}}}*/


; Phi, J, T

y4m2 = "201607"
path = "/data2/eraint/1.5x1.5/"

; Read file
; fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
; U  = fU->U_GDS0_ISBL 
; fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
; V  = fV->V_GDS0_ISBL 
fT = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T  = fT->T_GDS0_ISBL    ; Temperature [K]
fZ = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Z  = fZ->Z_GDS0_ISBL    ; Geopotential
fJ = addfile(path+"PHYpl/"+"Tphy.pl."+y4m2+".nc", "r")
J  = fJ->Tphy/(3600.*6) ; =J/Cp [K/6hour] -> [K/s]

lon = T&g0_lon_3
lat = T&g0_lat_2 
prs = T&lv_ISBL1*100    ; [hPa] -> [Pa]
a   = 6.37122e06        ; [m]
R   = 287.   ; 
pi  = 4*atan(1.)
d2r = pi/180.           ; Degree to radian
OMG = 7.2921e-5
P   = conform(T, prs, 1)

; Calculate Static stability 
; make only zonal average for practical reason
Tzm = dim_avg(T)
sig = static_stability(prs, T, 1, 0)
sig = (R/conform(sig, prs, 1))*sig  ; sig(t,p,y)

; Calculate f
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat
f({-3:3}) = f({3})          ; practical reasona
f0 = 1.e-4                  ; constant corioli parameter at 45N
; beta = beta_dfdy_rossby(lat,False)
; f0 = f-beta*(lat(2)-lat(1))*a                  ; constant corioli parameter
f2 = conform(sig, f, 2)

; Calculate geostrophic wind Ug, Vg
; Use f instead f0.
Ug = -(1./f2)*center_finite_diff_n(Z, a*lat*d2r, False, 1, 2)
Vg =  (1./f2)*center_finite_diff_n(Z, a*lon*d2r, False, 1, 3)
Vg = Vg/conform(Vg, cos(lat*d2r), 2)

; Calculate geostrophic vorticity
vo = uv2vr_cfd(Ug, Vg, lat, lon, 3)
; uv2vrf(Ug, Vg, vo)

; Forcing 1 : Absolute vorticity advection
cyclic = True ; grid is cyclic in longitude
longName = "Absolute vorticity advection"
units = "s**-2"
; abvo = -f0*advect_variable(u, v, x, gridType, longName, units, opt)
F1 = -f0*advect_variable_cfd( Ug,  Vg, conform(vo, f, 2)+vo, \
                             lat, lon, cyclic, longName, units, 0)
printVarSummary(F1)
exit

; Forcing 2 : Thickness advection
F2 = (-f0**2./sig)*advect_variable_cfd( Ug,  Vg,  R*T/P, \
                                       lat, lon, cyclic, longName, units, 0)
; Finite difference in log-P (d/dp -> 1/p*d/dlogp.)
logp = log(prs)
F2   = center_finite_diff_n(F2, logp, False, 1, 1)
F2   = F2/conform(F2, prs, 1)          ; take it back to pressure
printVarSummary(F2)
exit

; Forcing 3 : Heating
F3 = 

; Computes using lapsF function in ncl.
; For test:
;     ts    = shagC(tg)
;     n     = ispan(0,nlat-1,1)
;     nn1   = conform(ts,n*(n+1),1)
;     lap2s = -ts*nn1/(R^2)
;     lap2  = shsgC(lap2s,nlon)
; delZ = lapsF(Z)
















