; Calculate advection using spherical harmonics /*{{{*/
undef("advect_xy")
function advect_xy(u:numeric, v:numeric, x:numeric,gridType[1]:integer)
;
; compute advection of a variable: UV.GRADIENT(X) ; U*(dX/dlon) + V*(dX/dlat)
; (from advect_var in "$NCARG_ROOT/lib/ncarg/nclscripts/csm/contributed.ncl")
;
; Requires:
;    (1) grid be global cuz spherical harmonics are used
;    (2) required dimension order ([time,[lev,]lat,lon)
;    (3) the input grids *must* be ordered South==>North
; 
; Nomenclature:
;   u, v      - zonal and meridional wind components [m/s]
;               rightmost dimensions must be ([...,]lat,lon)
;   x         - scalar quantity to be advected 
;               rightmost dimensions must be ([...,]lat,lon)
;               eg: T, Z, divergence, vorticity, latent energy,..., whatever
;   gridType  - type of grid: 0=gaussian, 1=fixed
;
; Usage:
;   f = addfile ("foo.nc", "r")
;   u = f->U    ; (time,lev,lat,lon)
;   v = f->V
;   T = f->T
;               ; linear advection of temperature
;   Tadv = advect_xy(u,v,T,0)
;
local dimu, dimv, dimx, ranku, rankv, rankx, ier \
    , x_grad_lon, x_grad_lat, advect
begin
                               ; ERROR CHECKING
    dimu  = dimsizes(u)
    dimv  = dimsizes(v)
    dimx  = dimsizes(x)
    ranku = dimsizes(dimu)
    rankv = dimsizes(dimv)
    rankx = dimsizes(dimx)

    ier   = 0
    if (.not.(gridType.eq.0 .or. gridType.eq.1)) then
        print("advect_var: unrecognized gridType: only 0 and 1 allowed")
        print("            gridType="+gridType)
        ier = 1
    end if
    if (.not.(ranku.eq.rankv .and. ranku.eq.rankx)) then
        print("advect_var: all input arguments must have the same rank")
        print("            ranku="+ranku)
        print("            rankv="+rankv)
        print("            rankx="+rankx)
        ier = ier + 10
    end if
    if (.not.(all(dimu.eq.dimv) .and. all(dimu.eq.dimx))) then
        print("advect_var: all input arguments must have the same dimension sizes")
        print("            ranku="+ranku)
        print("            rankv="+rankv)
        print("            rankx="+rankx)
        ier = ier + 100
    end if
    if (rankx.lt.2) then
        print("advect_var: variable rank must be at least 2D: ([time,[lev,]lat,plon)")
        print("advect_var: rank="+rankx+"D")
        ier = ier + 1000
    end if

    dnam  = getvardims(x)             ; dimension names
    ynam  = dnam(rankx-2)             ; latitude dimension name

    if (.not.ismissing(ynam) .and. iscoord(x,ynam)) then
        yord  = isMonotonic(x&$ynam$)
        if (yord.le.0) then
            print("advect_var: grid is not in South-to-North order.")
            print("Note! Degree is set North-to-South for calculate advection.")
            u = u(:,:,::-1,:)
            v = v(:,:,::-1,:)
            x = x(:,:,::-1,:)
            ; ier = ier + 10000
        end if
    else
        print("advect_var: Warning: No named dimensions. Can not check grid order.")
    end if

    if (ier.ne.0) then
        print("advect_var fatal error(s) encountered: ier="+ier)
        exit
    end if

    x_grad_lon = new( dimx, typeof(x), getFillValue(x) )
    x_grad_lat = new( dimx, typeof(x), getFillValue(x) )

    if (gridType.eq.0) then
        gradsg  (x, x_grad_lon, x_grad_lat)
    end if
    if (gridType.eq.1) then
        gradsf  (x, x_grad_lon, x_grad_lat)
    end if
                                          ; advection

    uxgrad_lon = u*x_grad_lon 
    vxgrad_lat = v*x_grad_lat
    advect     = uxgrad_lon + vxgrad_lat
                                          ; add meta data

    if (.not.ismissing(ynam) .and. iscoord(x,ynam)) then
        yord  = isMonotonic(x&$ynam$)
        if (yord.le.0) then
            x = x(:,:,::-1,:)
            advect = advect(:,:,::-1,:)
            ; ier = ier + 10000
        end if
    end if

    copy_VarCoords(x, advect)

    return(advect)

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
T  = fT->T_GDS0_ISBL(:2,:,:,:)    ; Temperature [K]
fZ = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Z  = fZ->Z_GDS0_ISBL(:2,:,:,:)    ; Geopotential
fJ = addfile(path+"PHYpl/"+"Tphy.pl."+y4m2+".nc", "r")
J_cp  = fJ->Tphy(:2,:,:,:)/(3600.*6) ; =J/Cp [K/6hour] -> [K/s]

tim = T&initial_time0_hours
lon = T&g0_lon_3
lat = T&g0_lat_2 
prs = T&lv_ISBL1*100    ; [hPa] -> [Pa]
a   = 6.37122e06        ; [m]
R   = 287.              ; 
pi  = 4*atan(1.)
d2r = pi/180.           ; Degree to radian
OMG = 7.2921e-5
P   = conform(T, prs, 1)

; Calculate Static stability 
; make only zonal average for practical reason
Tzm = dim_avg(T)
sig = static_stability(prs, Tzm, 1, 0)
sig = (R/conform(sig, prs, 1))*sig  ; sig(t,p,y)

; Calculate f
f   = 2*OMG*sin(lat*d2r)    ; Coriolis force
f!0 = "lat"
f&lat = lat
f({0}) = 1.e-5          ; practical reason
                        ; Assume! f_equator is small value.
f4   = conform(T, f, 2)
beta = 2*OMG*cos(lat*d2r)/a

; Calculate geostrophic wind Ug, Vg
; Use f instead f0.
Ug = -(1./f4)*center_finite_diff_n(Z, a*lat*d2r, False, 1, 2)
Vg =  (1./f4)*center_finite_diff_n(Z, a*lon*d2r, False, 1, 3)
Vg = Vg/conform(Vg, cos(lat*d2r), 2)

; Calculate geostrophic vorticity
; vo = uv2vr_cfd(Ug, Vg, lat, lon, 3)
vo = uv2vrF(Ug, Vg)

; Forcing 1 : Absolute vorticity advection
grid = 1                      ; 1 : fixed grid, 0 : gaussian grid
avo  = f4+vo   ; Absolute vorcity
copy_VarCoords(T, avo)
F1 = -f4*advect_xy(Ug, Vg, avo, grid)   ; function in contribute.ncl

; Forcing 2 : Thickness advection
dZdp = R*T/P  ; -dgeo/dp = R*T/P
copy_VarCoords(T, dZdp)
sig  = where(sig .le. 0, 1.e-7, sig)   ; Assume! remove negative sigma.
sig4 = conform_dims(dimsizes(T), sig, (/0,1,2/))
F2   = (-f4^2./sig4)*advect_xy(Ug, Vg, dZdp, grid)
; optional (diff p)
F2   = -center_finite_diff_n(F2, prs, False, 1, 1)
; log-P
; logp = log(prs)       ; Finite difference in log-P (d/dp -> 1/p*d/dlogp.)
; F2   = center_finite_diff_n(F2, logp, False, 1, 1)
; F2   = -F2/P          ; take it back to pressure

; Forcing 3 : Heating
F3   = ((R/P)*f4^2./sig4*J_cp) 
; optional (diff p)
F3   = -center_finite_diff_n(F3, prs, False, 1, 1)
; log-P
; F3   = center_finite_diff_n(F3, logp, False, 1, 1)
; F3   = -F3/conform(F3, prs, 1)          ; take it back to pressure

; =============== Need to work================
;================================================================
exit

; Tendency term X = (dZ/dt)
; Computes using lapsF function in ncl.
; For test:
;     ts    = shagC(tg)
;     n     = ispan(0,nlat-1,1)
;     nn1   = conform(ts,n*(n+1),1)
;     lap2s = -ts*nn1/(R^2)
;     lap2  = shsgC(lap2s,nlon)
; delZ = lapsF(Z)
X    = center_finite_diff_n(Z, tim, False, 1, 0)
Xxy  = lapsF(X)     ; d(dX/dx)/dx + d(dX/dy)/dy

; Solve equation using SOR(Successive Over Relexation) method

; ============== test ==============
copy_VarCoords(T, F1)
copy_VarCoords(T, F2)
copy_VarCoords(T, F3)
copy_VarCoords(T, X)
copy_VarCoords(T, Xxy)
ofile = "./geotend.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->F1  = F1
fc->F2  = F2
fc->F3  = F3
fc->X   = X
fc->Xxy = Xxy
exit
















