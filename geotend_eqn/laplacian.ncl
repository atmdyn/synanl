load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_code.ncl"
load "$NCARG_ROOT/lib/ncarg/nclscripts/csm/gsn_csm.ncl"


a = addfile("fms.1500.nc","r")
R = 6371.0e3
set_sphere_radius(R)  ; fms value

nlat = 64
nlon = 128
temp = a->temp
tg = temp(0,39,:,:)
; 1 laplacian using function
lap1 = lapsG(tg)

; 2 laplacian using concept
ts = shagC(tg)
n  = ispan(0,nlat-1,1)
nn1 = conform(ts,n*(n+1),1)
lap2s = -ts*nn1/(R^2)
lap2  = shsgC(lap2s,nlon)

w=gsn_open_wks("pdf","laplacian")
res = False
p=gsn_csm_contour(w,lap1,res)
p=gsn_csm_contour(w,lap2,res)
p=gsn_csm_contour(w,lap1-lap2,res)
