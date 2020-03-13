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

;-- Set variables
; Define dt instead of relaxation coefficient for convenience
dt = 1.e9    ; 
resid = 9.e-8  ; critical value
f0 = 1.e-4
sig = 2.5e-6   ; [m^2 / Pa^2 / s^2] holton 4th 150p
; sig = 5.e-6   ; [m^2 / Pa^2 / s^2] holton 4th 150p

; Coordinate
lev = ispan(200,1000,50)*100. ; Pa
lat = ispan(-90,  90, 1)    ; 181
lon = ispan(  0, 359, 1)    ; 360
nlev = dimsizes(lev)
nlon = dimsizes(lon)

; Wavelength
Lz = 20000.     ; [Pa]
Ly = 10.
Lx = 10. 

; Earth radius [m]
R = get_sphere_radius()    ; <- update (2020.01.29)


;-- Forcing Part
f1 = new((/nlev,dimsizes(lat),nlon/), "float")
hgt3 = conform(f1, lev, 0)
lat3 = conform(f1, lat, 1)
lon3 = conform(f1, lon, 2)

;=================== dev version =================== /*{{{*/
; Gaussian forcing 
A = 10./86400.   ; amplitude 10 [K/day]
; 1. One cell forcing
;  -1. Bottom forcing
; f1 =  A* exp(- ( ( lat3      /Ly )^2 \
;                + ((lon3-180) /Lx )^2 \    ; 90
;                + ((hgt3    ) /Lz )^2 ) )
;  -2. Z- Center forcing
f1 =  A* exp(- ( ( lat3      /Ly )^2 \
               + ((lon3-180) /Lx )^2 \    ; 90
               + ((hgt3-50000) /Lz )^2 ) )
;
; 2. Two cell forcing
;  -1. Equator 2 forcing (90: warming, 270: cooling)
; f1 =  A* exp(- ( ( lat3     /Ly )^2 \
;                + ((lon3-90) /Lx )^2 \    ; 90
;                + ((hgt3-10) /Lz )^2 ) )
; f1 = -A* exp(- ( ( lat3     /Ly )^2 \
;                + ((lon3-270)/Lx )^2 \    ; 270
;                + ((hgt3-10) /Lz )^2 ) ) + f1
;  -2. Polar 2 forcing 
; f1 =  A* exp(- ( ((lat3+90) /Ly )^2 \
;                + ((lon3-90) /Lx )^2 \    ; 90
;                + ((hgt3   ) /Lz )^2 ) )
; f1 = -A* exp(- ( ((lat3+90)    /Ly )^2 \
;                + ((lon3-270)/Lx )^2 \    ; 270
;                + ((hgt3   ) /Lz )^2 ) ) + f1
;  -3. SN 2 forcing 
; f1 =  A* exp(- ( ((lat3-90)  /Ly )^2 \
;                + ((hgt3   )  /Lz )^2 ) )
;                ; + ((lon3-180) /Lx )^2 \    ; 90
; f1 = -A* exp(- ( ((lat3+90)  /Ly )^2 \
;                + ((hgt3   )  /Lz )^2 ) ) + f1
;                ; + ((lon3-180) /Lx )^2 \    ; 270
;=================================================== /*}}}*/
b_lm = shaeC(f1)   ; [2] x [21] x [181] x [181]

;-- Solver Part
dp = lev(1:nlev-1) - lev(0:nlev-2)
sig := conform(dp, sig, 0)  ; index = k + 0.5

a_lm = todouble(b_lm*0.)    ; Because R is double
; f1 = f1*0.                ; <- Forcing test
; a_lm(0,0,1,1) = 1.        ; bottom boundary
n = ispan(0,dimsizes(lat)-1,1)
nn1 = conform(a_lm(:,0,:,:),n*(n+1),1)  ; [2] x [l] x [m]

err = 1.d5
begTime = get_cpu_time()
do while (err .gt. resid)
    dummy = a_lm
    do k = 1, dimsizes(lev)-2   ; omit boundary
        relax1 = (f0^2)*dt/(sig(k  )*dp(k  )*0.5*(lev(k+1)-lev(k-1)))
        relax2 = (f0^2)*dt/(sig(k-1)*dp(k-1)*0.5*(lev(k+1)-lev(k-1)))
        laps = (nn1/R^2)*dt*a_lm(:,k,:,:) ;*100 ; for test
        laps(:,0,0) = 0.        ; No mean value for laplacian
        a_lm(:,k,:,:) = a_lm(:,k,:,:)                          \
                        - laps                                 \
                        - (relax1+relax2)*a_lm(:,k,:,:)        \
                        + relax1*a_lm(:,k+1,:,:)               \
                        + relax2*a_lm(:,k-1,:,:)               \
                        + dt*b_lm(:,k,:,:)
    end do    
    err = avg(abs(a_lm - dummy))
    print(err+"")
end do    
solution = shseC(a_lm, nlon)
print("Relaxation time: " + (get_cpu_time() - begTime)+" [s]")

ofile = "./3d_output.nc"
system("rm -f "+ofile)
f1!0 = "hgt"
f1!1 = "lat"
f1!2 = "lon"
copy_VarCoords(f1, solution)
fc = addfile(ofile, "c")
fc->forcing  = f1
fc->solution = solution
