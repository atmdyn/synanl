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
; Relaxation coefficient : relax = (f0/N)^2*(2*dt/dz^2)
relax = 0.004    ; = 1 : normal relaxation
               ; > 1 : over relaxationa
resid = 1.e-4  ; critical value
f0 = 1.e-4
N  = 1.e-1

; Coordinate
lev = ispan(  0,  20, 1)    ; [km]
lat = ispan(-90,  90, 1)    ; 181
lon = ispan(  0, 359, 1)    ; 360
nlon = dimsizes(lon)

; Wavelength
Lz = 2.     ; [km]
Ly = 10.
Lx = 10. 

; Earth radius [m]
R = get_sphere_radius()    ; <- update (2020.01.29)


;-- Forcing Part
f1 = new((/21,dimsizes(lat),dimsizes(lon)/), "float")
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
; f1 =  A* exp(- ( ( lat3      /Ly )^2 \
;                + ((lon3-180) /Lx )^2 \    ; 90
;                + ((hgt3- 10) /Lz )^2 ) )
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
f1 =  A* exp(- ( ((lat3-90)  /Ly )^2 \
               + ((hgt3   )  /Lz )^2 ) )
               ; + ((lon3-180) /Lx )^2 \    ; 90
f1 = -A* exp(- ( ((lat3+90)  /Ly )^2 \
               + ((hgt3   )  /Lz )^2 ) ) + f1
               ; + ((lon3-180) /Lx )^2 \    ; 270
;=================================================== /*}}}*/
b_lm = shaeC(f1)   ; [2] x [21] x [181] x [181]

;-- Solver Part
dz = 1000.  ; [m]
dt = ((N/f0)^2*(dz)^2)/2.*relax

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
