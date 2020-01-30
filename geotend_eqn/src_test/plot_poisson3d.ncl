a = addfile("./3d_output.nc", "r")
f = a->forcing
s = a->solution

; index
z = 20
y = 90
x = 180

wks = gsn_open_wks("x11", "./fig/nonamed")
res = True
res@cnFillOn = True
; Forcing
res@gsnLeftString = "Forcing x-y"
plt = gsn_csm_contour(wks, f(z,:,:), res)

res@gsnLeftString = "Forcing z-y"
plt = gsn_csm_contour(wks, f(:,:,x), res)

res@gsnLeftString = "Forcing z-x"
plt = gsn_csm_contour(wks, f(:,y,:), res)

; Solution
res@gsnLeftString = "Solution x-y"
plt = gsn_csm_contour(wks, s(z,:,:), res)

res@gsnLeftString = "Solution z-y"
plt = gsn_csm_contour(wks, s(:,:,x), res)

res@gsnLeftString = "Solution z-x"
plt = gsn_csm_contour(wks, s(:,y,:), res)

; system("crop ./fig/nonamed.pdf")
; system("montage ./fig/*.png -tile 3x2 -geometry 1000 out.pdf")
