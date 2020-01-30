; ------------------------------------------------ ;
; Created on  Thu Jan 30  11:01:24 2020
; @author : Joowan Kim, Hyeong-Gyu Kim. Hyeon-Seok Do
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
;  - [ ( laps )T = -F ]
;  - laps : laplacian 2d,  T : Temperature, 
;    F    : Forcing
;  - Laplacian part solve using spherical harmonics.
;       -l(l+1)
;    ΣΣ _______ * a_lm * Y_lm = ΣΣ b_lm * Y_lm
;         R^2
;   (See google drive 
;    /project_doc/TODO_Changma/etc/Spheical Harmonics.pdf)
;
; 4. Update
;  - 2020.01.29 : ERROR
;                toint(), tofloat() make segmentation fault
;                R = tofloat(get_sphere_radius())
; ------------------------------------------------ ;

; Coordinate
lat = ispan(-90, 90, 1)
lon = ispan(0, 359, 1)

Lx = 10.
Ly = 10.
R = get_sphere_radius()

f1 = new((/181,360/), "float") ; forcing
lat2 = conform(f1, lat, 0)
lon2 = conform(f1, lon, 1)

; Gaussian forcing (90: warming, 270: cooling)
f1 = exp(-((lat2/Ly)^2 + ((lon2-90)/Lx)^2))       ; 90
f1 = f1 - exp(-((lat2/Ly)^2 + ((lon2-270)/Lx)^2)) ; 270
b_lm = shaeC(-f1)   ; [2] x [181] x [181]

a_lm = todouble(b_lm)*0.
do l = 1, dimsizes(lat)-1
    a_lm(:,l,:) = -(R^2)/(l*(l+1))*b_lm(:,l,:)
end do    

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
plt = gsn_csm_contour(wks, f1, res)
plt = gsn_csm_contour(wks, a_lm(0,:,:), res)
plt = gsn_csm_contour(wks, shseC(a_lm, dimsizes(lon)), res)
