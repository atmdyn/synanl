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
lev1 = 200  ; [hPa]
lev2 = 500  ;
cp   = 1005.  ; [J/kg/K]
casetime = "2016070500"

ct = stringtochar(casetime)
y4 = chartostring(ct(0:3))
m2 = chartostring(ct(4:5))
d2 = chartostring(ct(6:7))
h2 = chartostring(ct(8:9))
y4m2 = y4+m2

; Read data
f = addfile("apparent_heat_moist.nc", "r")
Q1 = f->Q1
Q2 = f->Q2
t = Q1&initial_time0_hours
unit = t@units
ctim = cd_inv_calendar(toint(y4), toint(m2), toint(d2), toint(h2), 0, 0, unit, 0)
tim1 = cd_inv_calendar(toint(y4), toint(m2), toint(d2)-1, toint(h2), 0, 0, unit, 0)
tim2 = cd_inv_calendar(toint(y4), toint(m2), toint(d2)+1, toint(h2), 0, 0, unit, 0)

; Calculate for budget anal
data1 = Q1({tim1:tim2},:,:,:) / cp * 86400.
data1 := dim_avg_n(data1, 0)    ; [K/day]
copy_VarCoords(Q1({ctim},:,:,:), data1)
data2 = Q2({tim1:tim2},:,:,:) / cp * 86400.
data2 := dim_avg_n(data2, 0)    ; [K/day]
copy_VarCoords(Q2({ctim},:,:,:), data2)

; For check
ftphy = addfile("/data2/eraint/1.5x1.5/PHYpl/Tphy.pl."+y4m2+".nc", "r")
Tphy = ftphy->Tphy
data3 = Tphy({tim1:tim2},:,:,:) * 4.
data3 := dim_avg_n(data3, 0)    ; [K/day]
copy_VarCoords(Q1({ctim},:,:,:), data3)

wks = gsn_open_wks("x11", "./apparent_"+casetime)
res = True
res@cnFillOn = True
res@gsnAddCyclic = False

res@mpMinLatF = 20
res@mpMaxLatF = 50
res@mpMinLonF = 100
res@mpMaxLonF = 150
;res@cnLevelSelectionMode = "ExplicitLevels"
;res@cnLevels = ispan(-30,30,1)

res@gsnLeftString = "Q1("+lev1+"hPa)"
res@gsnRightString = "K/day"
plt = gsn_csm_contour_map(wks, data1({lev1},:,:), res)
res@gsnLeftString = "Q2("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, data2({lev1},:,:), res)
res@gsnLeftString = "Q("+lev1+"hPa)"
plt = gsn_csm_contour_map(wks, data3({lev1},:,:), res)

res@gsnLeftString = "Q1("+lev2+"hPa)"
plt = gsn_csm_contour_map(wks, data1({lev2},:,:), res)
res@gsnLeftString = "Q2("+lev2+"hPa)"
plt = gsn_csm_contour_map(wks, data2({lev2},:,:), res)
res@gsnLeftString = "Q("+lev2+"hPa)"
plt = gsn_csm_contour_map(wks, data3({lev2},:,:), res)


;======================================================


; average
Q1_prof = dim_avg(dim_avg(data1(:,{34:42},{125:130})))
Q2_prof = dim_avg(dim_avg(data2(:,{34:42},{125:130})))
Q_prof  = dim_avg(dim_avg(data3(:,{34:42},{125:130})))

lres = True
lres@trYReverse = True
lres@trYLog = True
lres@trYMinF = 1
lres@xyLineColors = (/"red", "blue", "black"/)
lres@xyDashPatterns = (/0,0,0/)
lres@xyLineThicknessF = 3.
lres@tmYLMode = "Explicit"
lres@tmYLValues = (/ 1000, 850, 500, 300, 100, 30, 10, 3, 1 /)
lres@tmYLLabels = lres@tmYLValues

lres@gsnLeftString = "Q1(Red), Q2(Blue), Q(Black)"
res@gsnRightString = "K/day"

prs = Q1&lv_ISBL1
plt = gsn_csm_xy(wks, (/Q1_prof,Q2_prof,Q_prof/), prs, lres)





