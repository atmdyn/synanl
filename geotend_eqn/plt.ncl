f = addfile("./output.nc", "r")
s = f->solution

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
plt = gsn_csm_contour(wks, s(0,:,:), res)
plt = gsn_csm_contour(wks, s(:,:,90), res)
plt = gsn_csm_contour(wks, s(:,:,270), res)
