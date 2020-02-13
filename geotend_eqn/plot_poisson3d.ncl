a = addfile("./geotend.nc", "r")
f = a->forcing
s = a->solution
lat = ispan(-90,  90, 1)    ; 181
lat@units = "degree_N"
lon = ispan(  0, 359, 1) 
lon@units = "degree_E"
f&lat = lat
f&lon = lon
copy_VarCoords(f, s)

; index
z = 5
y = 90
x = 180

wks = gsn_open_wks("x11", "./fig/nonamed")
res = True
res@cnFillOn = True
; Forcing
res@gsnLeftString = "Forcing x-y"
; plt = gsn_csm_contour(wks, f(z,:,:), res)
res@mpCenterLonF = 180
plt = gsn_csm_contour_map(wks, f(z,:,:), res)
; plt = gsn_csm_contour_map_polar(wks, f(z,:,:), res)

res@gsnLeftString = "Forcing z-y"
plt = gsn_csm_contour(wks, f(:,:,x), res)

res@gsnLeftString = "Forcing z-x"
plt = gsn_csm_contour(wks, f(:,y,:), res)

; Solution
res@gsnLeftString = "Solution x-y"
; plt = gsn_csm_contour(wks, s(z,:,:), res)
plt = gsn_csm_contour_map(wks, s(z,:,:), res)
; res@gsnPolar = "SH"
; plt = gsn_csm_contour_map_polar(wks, s(z,:,:), res)

res@gsnLeftString = "Solution z-y"
plt = gsn_csm_contour(wks, s(:,:,x), res)

res@gsnLeftString = "Solution z-x"
plt = gsn_csm_contour(wks, s(:,y,:), res)

; system("crop ./fig/nonamed.pdf")
; system("montage ./fig/*.png -tile 3x2 -geometry 1000 out.pdf")
