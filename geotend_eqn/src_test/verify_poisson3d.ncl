; ------------------------------------------------ ;
; Created on  Thu Jan 30  11:01:24 2020
; @author : Joowan Kim, Hyeong-Gyu Kim, Hyeon-Seok Do
;
; 3D Poisson's equation solver (test version)
; 1. Filename : poisson3d.ncl
;
; 2. Variable
;  - Input  : boundary condition
;  - Output : output
;
; 3. Description
;  - Compute poisson's equation 
;  - [ ( laps + (f0/N)^2*(d2/dz2) )T = -F ]
;  - laps : laplacian 2d,   f0 : Corioli at 45N
;    T    : Temperature,    N  : Bouyancy frequency
;    F    : Forcing
;  - Laplacian part solve using spherical harmonics.
;   (See google drive 
;    /project_doc/TODO_Changma/etc/Spheical Harmonics.pdf)
;  - Vertical diffusion part is solved using Relaxaion method.
;
; 4. Update
;  - 2020.01.29 : ERROR
;                toint(), tofloat() make segmentation fault
;                R = tofloat(get_sphere_radius())
; ------------------------------------------------ ;

load "poisson_solver.ncl"

;-- Set variables
; Define dt instead of relaxation coefficient for convenience
dt = 1.e9    ; 
resid = 1.e-13  ; critical value
f0 = 1.e-4
sig = 2.5e-6   ; [m^2 / Pa^2 / s^2] holton 4th 150p

; Coordinate
lev = ispan(200,1000,50)*100. ; Pa
lat = ispan(-90,  90, 1)    ; 181
lon = ispan(  0, 359, 1)    ; 360
nlev = dimsizes(lev)
nlon = dimsizes(lon)

; Wavelength
Lz = 40000.     ; [Pa]
Ly = 40.
Lx = 40. 

; Earth radius [m]
R = get_sphere_radius()    ; <- update (2020.01.29)


;-- Solution Part
f1 = new((/nlev,dimsizes(lat),nlon/), "float")
hgt3 = conform(f1, lev, 0)
lat3 = conform(f1, lat, 1)
lon3 = conform(f1, lon, 2)

;=================== dev version =================== /*{{{*/
; Gaussian forcing 
A  = 200./86400.   ; amplitude 200 [m/day]
S1 =  A* exp(- ( ( lat3      /Ly )^2 \
               + ((lon3-180  ) /Lx )^2 \
               + ((hgt3-60000) /Lz )^2 ) )
;=================================================== /*}}}*/

; compute laplacian
laph = lapsF(S1)
dp   = 5000.
iz   = ispan(1, nlev-2, 1)
lapv = S1
lapv(iz,:,:) = ( (f0^2)/sig)*(S1(iz+1,:,:) \
                          - 2*S1(iz  ,:,:) \
                            + S1(iz-1,:,:) )/(dp*dp)
                           
lapv(0,:,:) = 0.
lapv(nlev-1,:,:) = 0.
lap  = laph + lapv ; <- forcing
lap!0 = "lev"
lap!1 = "lat"
lap!2 = "lon"
lap&lev = lev
lap&lat = lat
lap&lon = lon

; Inverse laplacian using poisson solver
opt = True
opt@lev_units = "Pa"
opt@f0 = f0
opt@sig = sig
opt@thres = resid
opt@dt = dt
opt@upper_bound = S1(0,:,:)
opt@lower_bound = S1(nlev-1,:,:)
sol = poisson_solver(lap, opt)

copy_VarCoords(lap, S1)

ofile = "./3d_output.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->X1 = S1
fc->X2 = sol
fc->F  = lap
