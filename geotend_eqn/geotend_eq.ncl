; /*{{{*/
undef("laplcaian")
function laplcaian(Y)
; argv = [time] x [lev] x [lat] x [lon]
begin
    
    dx  = Y&lon
    Yxx = (Y(:,:,:,ix+1) - 2*Y(:,:,:,ix) + Y(:,:,:,ix-1)) / (dx**2)

    return 0
end ; /*}}}*/


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

; Calculate Static stability 
; make only zonal average for practical reason
Tzm = dim_avg(T)
sig = static_stability(prs, T, 1, 0)
sig = (R/conform(sig, prs, 1))*sig  ; sig(t,p,y)

; Calculate f
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat
f({-3:3}) = f({3})          ; practical reason

; Calculate geostrophic wind Ug, Vg
; Use f instead f0.
Ug = -(1./f)*center_finite_diff_n(Z, a*lat*d2r, False, 1, 2)
Vg =  (1./f)*center_finite_diff_n(Z, a*lon*d2r, False, 1, 3)
Vg = Vg/conform(Vg, cos(lat*d2r), 2)
printVarSummary(Vg)

; Calculate geostrophic vorticity
















