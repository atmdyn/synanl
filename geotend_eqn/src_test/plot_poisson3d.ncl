a = addfile("./3d_output.nc", "r")
s1 = a->X1    ; maded
F  = a->F     ; calculated
s2 = a->X2    ; calculated
lat = ispan(-90,  90, 1)    ; 181
lat@units = "degree_N"
lon = ispan(  0, 359, 1) 
lon@units = "degree_E"
F&lat = lat
F&lon = lon
copy_VarCoords(F, s1)
copy_VarCoords(F, s2)

diff = s1 - s2
copy_VarCoords(s1, diff)

; index
z = 5
y = 90
x = 180

; 1. scale (Lx, Lz, ..)
; 2. near surface
; 3. multiple synoptic pattern

wks = gsn_open_wks("x11", "./fig/nonamed")
res = True
res@cnFillOn = True
res1 = res
res1@mpCenterLonF = 180
res2 = res
res2@trYReverse = True
; res2@trYLog = True

;-- S1
res1@gsnLeftString = "X1 x-y"
plt = gsn_csm_contour_map(wks, diff(z,:,:), res1)

res2@gsnLeftString = "X1 z-y"
plt = gsn_csm_contour(wks, diff(:,:,x), res2)

res2@gsnLeftString = "X1 z-x"
plt = gsn_csm_contour(wks, diff(:,y,:), res2)
exit
;-- S1
res1@gsnLeftString = "X1 x-y"
plt = gsn_csm_contour_map(wks, s1(z,:,:), res1)

res2@gsnLeftString = "X1 z-y"
plt = gsn_csm_contour(wks, s1(:,:,x), res2)

res2@gsnLeftString = "X1 z-x"
plt = gsn_csm_contour(wks, s1(:,y,:), res2)

;-- Forcing
res1@gsnLeftString = "F x-y"
plt = gsn_csm_contour_map(wks, F(z,:,:), res1)

res2@gsnLeftString = "F z-y"
plt = gsn_csm_contour(wks, F(:,:,x), res2)

res2@gsnLeftString = "F z-x"
plt = gsn_csm_contour(wks, F(:,y,:), res2)

;-- S2
res1@gsnLeftString = "X2 x-y"
plt = gsn_csm_contour_map(wks, s2(z,:,:), res1)

res2@gsnLeftString = "X2 z-y"
plt = gsn_csm_contour(wks, s2(:,:,x), res2)

res2@gsnLeftString = "X2 z-x"
plt = gsn_csm_contour(wks, s2(:,y,:), res2)

