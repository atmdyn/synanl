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
y4m2 = "201601"


; Read data
f1 = addfile("geotend.nc", "r")
; X  = f1->X
; sol = f1->solution

path = "/data2/eraint/1.5x1.5/"
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
PHI  = fin->Z_GDS0_ISBL(0:2,{200:1000},:,:)    ; Geopotential [m**2 s**-2]
X    = (PHI(2,:,:,:) - PHI(0,:,:,:)) / 21600. /2.
; X    = center_finite_diff_n(PHI(::4,:,:,:), 21600., False, 1, 0)
; X   = dim_avg_n_Wrap(f1->X (:,:,:,:),0)
sol = dim_avg_n_Wrap(f1->solution(1:1,:,:,:),0)
; X := X(1,:,:,:)
copy_VarCoords(sol, X)


wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
res@cnLinesOn = False
res@cnFillPalette = "amwg256"
res@cnLevelSelectionMode = "ExplicitLevels"
res@cnLevels = ispan(-40,40,2)*1.e-3
; res@gsnAddCyclic = False

; res@mpMinLatF = 20
; res@mpMaxLatF = 50
; res@mpMinLonF = 100
; res@mpMaxLonF = 150

; res@gsnRightString = "K/day"
res@gsnLeftString = "X ("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, X({lev1},:,:), res)

; res@cnLevels := fspan(-4,4,10)*1.e-13
; res@gsnLeftString = "F1("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F1({lev1},:,:), res)
;
; res@gsnLeftString = "F2("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F2({lev1},:,:), res)
;
; res@cnLevels := fspan(-4,4,10)*1.e-15
; res@gsnLeftString = "F3("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F3({lev1},:,:), res)

; res@cnLevels := ispan(-24,24,4)*1.e-3
res@gsnLeftString = "solution ("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, sol({lev1},:,:), res)
