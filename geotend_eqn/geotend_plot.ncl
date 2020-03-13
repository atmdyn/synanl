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
; f = addfile("geotend.nc", "r")
flist = systemfunc("ls -1 geotend_*.nc")
f = addfiles(flist, "r")
ListSetType(f, "join")
; F1 = dim_avg_n_Wrap(f[:]->F1,0)
; F2 = dim_avg_n_Wrap(f[:]->F2,0)
; F3 = dim_avg_n_Wrap(f[:]->F3,0)
X  = dim_avg_n_Wrap(f[:]->X ,0)
sol = dim_avg_n_Wrap(f[:]->solution,0)


wks = gsn_open_wks("pdf", "nonamed")
res = True

res@cnFillOn = True
res@cnLevelSelectionMode = "ExplicitLevels"
; res@gsnAddCyclic = False

; res@mpMinLatF = 20
; res@mpMaxLatF = 50
; res@mpMinLonF = 100
; res@mpMaxLonF = 150
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = ispan(-30,30,1)

; res@gsnRightString = "K/day"
res@gsnLeftString = "X ("+lev1+"hPa)"
y4m2 = "201602"
path = "/data2/eraint/1.5x1.5/"
time_idx = 4
fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
Zg   = fin->Z_GDS0_ISBL(0:11,{200:1000},:,:)    ; Geopotential
dims = dimsizes(Zg)
Zd   = dim_avg_n_Wrap(reshape(Zg, (/3,4,dims(1),dims(2),dims(3)/)),1)
X   := center_finite_diff_n(Zd, 86400., False, 1, 0)
copy_VarCoords(Zg(0,:,:,:), X(0,:,:,:))
plt = gsn_csm_contour_map(wks, X(1,{lev1},:,:), res)
res@cnLevels = ispan(-24,24,4)*1.e-3

; res@gsnLeftString = "F1("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F1({lev1},:,:), res)
;
; res@gsnLeftString = "F2("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F2({lev1},:,:), res)
;
; res@gsnLeftString = "F3("+lev1+"hPa)"
; plt = gsn_csm_contour_map(wks, F3({lev1},:,:), res)

res@cnLevels := ispan(-6,6,1)*1.e-2
res@gsnLeftString = "solution("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, sol({lev1},:,:), res)
