; ------------------------------------------------ ;
; Created on  Thu Jan 02  18:37:01 2020
; @author : Joowan Kim, Hyeong-Gyu Kim. Hyeon-Seok Do
;
;-- Description
;  - Check about apparent result(./apparent_calc.ncl).
;
; 
;-- Update
;  - 2020-01-02 : Maked by ADLab
; ------------------------------------------------ ;

; Set variable
lev1 = 500  ;
cp   = 1005.  ; [J/kg/K]
casetime = "2016070500"

ct = stringtochar(casetime)
y4 = chartostring(ct(0:3))
m2 = chartostring(ct(4:5))
d2 = chartostring(ct(6:7))
h2 = chartostring(ct(8:9))
y4m2 = y4+m2

; Read data
f = addfile("geotend.nc", "r")
F1 = f->F1 (:,:,:)
F2 = f->F2 (:,:,:)
F3 = f->F3 (:,:,:)
X  = f->X  (:,:,:)
sol  = f->solution  (:,:,:)


wks = gsn_open_wks("x11", "nonamed")
res = True

res@cnFillOn = True
res@gsnAddCyclic = False

; res@mpMinLatF = 20
; res@mpMaxLatF = 50
; res@mpMinLonF = 100
; res@mpMaxLonF = 150
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = ispan(-30,30,1)

res@gsnLeftString = "X ("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, X({lev1},:,:), res)
res@gsnLeftString = "F1("+lev1+"hPa)"
; res@gsnRightString = "K/day"
plt = gsn_csm_contour_map(wks, F1({lev1},:,:), res)
res@gsnLeftString = "F2("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, F2({lev1},:,:), res)
res@gsnLeftString = "F3("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, F3({lev1},:,:), res)
res@gsnLeftString = "solution("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, sol({lev1},:,:), res)
