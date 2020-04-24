load "/home/gue/Dropbox/work/TODO/synopanal/geotend_eqn/src_test/poisson_solver.ncl"
kc1 = asciiread("/home/hyey/coldsurge/2019/synop/first_cold/clust/runavg/new_cl/kc1.txt",28,"integer")
fi = addfile("/work1/hyey/data/ecmwf/streamf/qgpv/v2/streamf_cold.b1.nc","r")
sf = dim_avg_n_Wrap(fi->sf,0)
id = asciiread("/home/hyey/coldsurge/2019/synop/first_cold/clust/runavg/new_cl/kc2.txt",33,"integer")
dir  = "/work1/hyey/data/ecmwf/"
tim1 = ispan(1,38*12,1)-1
tim1@units = "monthes since 1981-01-01 0:0:0"
case = asciiread("/home/hyey/coldsurge/2018/quant/getcase37.txt",(/61,2/),"float")
event = toint(case(:,0))
event@units = "days since 1981-01-01 0:0:0"
flist1 = dir+"Tpl/Tpl."+cd_calendar(tim1,-1)+".day.nc"
fhead1 = addfiles(flist1,"r")
T = fhead1[:]->temp(event,:,:,:)
T2 = T(kc1,:,:,:)

fi = addfile("/work1/hyey/data/ecmwf/streamf/qgpv/v2/tqgpv_kc1.85.b1.nc","r")
mqgpv= fi->qgpv
sh = fi->SH
cy = fi->CY
sh_str = fi->SH_str
cy_str = fi->CY_str
str = fi->str
lat = fi->lat
lon = fi->lon
lev = fi->level
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)

input1 = sh
input2 = cy
input3 = mqgpv
;input = mqgpv
;printVarSummary(input)


;sf_z = center_finite_diff_n(sf,lev*100,False,0,0)
;printMinMax( sf_z,False)

sf_b = (sf(nlev-1,:,:)-sf(nlev-2,:,:))/((lev(nlev-1)-lev(nlev-2))*100.)*10
sf_t = (sf(1,:,:)-sf(0,:,:))/((lev(1)-lev(0))*100.)
copy_VarCoords(sf(0,:,:),sf_t)
copy_VarCoords(sf(0,:,:),sf_b)

printMinMax( sf_t,False)
printMinMax( sf_b,False)

;w = gsn_open_wks("x11","ee")
r = True
r@cnFillOn = True
;p = gsn_csm_contour_map(w,sf_b,r)

opt = True
opt@lev_units = "hPa"
opt@bound_type = "D"
opt@T0 = dim_avg_n_Wrap(dim_avg_Wrap(dim_avg_Wrap(T(:,:,:,:))),0)
opt@dt = 6.e8
; opt@thres = 0.9
; opt@f0 = 1.e-5
;opt@upper_bound = sh_str(0,:,:)
;opt@lower_bound = sh_str(nlev-1,:,:)
;result1 = poisson_solver(input1, opt)

;opt@upper_bound = cy_str(0,:,:)
;opt@lower_bound = cy_str(nlev-1,:,:)
;result2 = poisson_solver(input2, opt)
opt@upper_bound = sf_t;sf(0,:,:)
opt@lower_bound = sf_b;sf(nlev-1,:,:)
result3 = poisson_solver(input3, opt)

system("rm -f ./test_n10.nc")
fc = addfile("./test_n10.nc", "c")
;ystem("rm -f ./output/output_kc1.85.b1.allbd.nc")
;'fc = addfile("./output/output_kc1.85.b1.allbd.nc", "c")
;fc->sh_out = result1
;fc->cy_out = result2
fc->sf_out = result3
;fc->sf   = sf
