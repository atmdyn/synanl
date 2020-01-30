; ------------------------------------------------ ;
; Created on  Thu Jan 30  11:01:24 2020
; @author : Hyeong-Gyu Kim.
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

;-- Set variables
; Relaxation coefficient : relax = (f0/N)^2*(2*dt/dz^2)
relax = 0.2    ; = 1 : normal relaxation
               ; > 1 : over relaxationa
resid = 1.e-8  ; critical value
f0 = 1.e-4
N  = 1.e-2

; Coordinate
lev = ispan(  0,  20, 1)   ; [km]
lat = ispan(-90,  90, 1)
lon = ispan(  0, 359, 1)
nlon = dimsizes(lon)

; Wavelength
Lz = 2.     ; [km]
Ly = 10.
Lx = 10. 

; Earth radius [m]
R = get_sphere_radius()    ; <- update (2020.01.29)


;-- Forcing Part
f1 = new((/21,181,360/), "float")
hgt3 = conform(f1, lev, 0)
lat3 = conform(f1, lat, 1)
lon3 = conform(f1, lon, 2)

; Gaussian forcing (90: warming, 270: cooling)
A = 10./86400.   ; amplitude 10 [K/day]
f1 =  A* exp(- ( ( lat3     /Ly )^2 \
               + ((lon3-90) /Lx )^2 \    ; 90
               + ((hgt3-10) /Lz )^2 ) )
f1 = -A* exp(- ( ( lat3     /Ly )^2 \
               + ((lon3-270)/Lx )^2 \    ; 270
               + ((hgt3-10) /Lz )^2 ) ) + f1
f1 = f1*0.         ; <- for forcing test
b_lm = shaeC(f1)   ; [2] x [21] x [181] x [181]


;-- Solver Part
dz = 1000.  ; [m]
dt = ((N/f0)^2*(dz)^2)/2.*relax

a_lm = todouble(b_lm*0.)    ; Because R is double
a_lm(0,0,2,2) = 1.          ; bottom boundary
n = ispan(0,dimsizes(lat)-1,1)
nn1 = conform(a_lm(:,0,:,:),n*(n+1),1)  ; [2] x [l] x [m]

err = 1.d5
do while (err .gt. resid)
    dummy = a_lm
    do k = 1, dimsizes(lev)-2   ; omit boundary
        laps = (nn1/R^2)*dt*a_lm(:,k,:,:) 
        laps(:,0,0) = 0.        ; No mean value for laplacian
        a_lm(:,k,:,:) = a_lm(:,k,:,:) - relax*a_lm(:,k,:,:)            \
                        - laps                                         \
                        + relax/2.*(a_lm(:,k+1,:,:) + a_lm(:,k-1,:,:)) \
                        + dt*b_lm(:,k,:,:)
    end do    
    err = avg(abs(a_lm - dummy))
    print(err+"")
end do    

solution = shseC(a_lm, nlon)

ofile = "./3d_output.nc"
system("rm -f "+ofile)
f1!0 = "hgt"
f1!1 = "lat"
f1!2 = "lon"
copy_VarCoords(f1, solution)
fc = addfile(ofile, "c")
fc->forcing  = f1
fc->solution = solution
