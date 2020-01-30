; ------------------------------------------------ ;
; Created on  Thu Jan 30  10:57:36 2020
; @author : Joowan Kim, Hyeong-Gyu Kim, Hyeon-Seok Do
;
; 1. Filename : verify.ncl
;
; 2. Variable
;  - Input  : result of "poisson3d.ncl"
;             = (laplacian + (f0/N)^2*(d2/dz2))T
;  - Output : -F (Forcing)
;
; 3. Description
;  - Compute forcing for verify
; ------------------------------------------------ ;

f = addfile("./output.nc", "r")
s = f->solution

f0 = 1.e-4
N  = 1.e-2
dz = 1000 ; m
Ts = lapsF(s)

Tz = s*0.
do k = 1, 19
    Tz(k,:,:) = (f0/N)^2*(s(k+1,:,:)+s(k-1,:,:)-2*s(k,:,:))/(dz*dz)
end do    
f1 = -(Ts+Tz)

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
plt = gsn_csm_contour(wks, f1(10,:,:), res)

