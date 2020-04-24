id3 = asciiread("/home/hyey/coldsurge/2019/synop/first_cold/clust/runavg/new_cl/kc1.txt",28,"integer")

case = asciiread("/home/hyey/coldsurge/2018/quant/getcase37.txt",(/61,2/),"float")
event = toint(case(:,0))-1
event@units = "days since 1981-01-01 0:0:0"

dir  = "/work1/hyey/data/ecmwf/"
tim1 = ispan(1,38*12,1)-1
tim1@units = "monthes since 1981-01-01 0:0:0"
flist1 = dir+"Tpl/Tpl."+cd_calendar(tim1,-1)+".day.nc"
fhead1 = addfiles(flist1,"r")
tpl = fhead1[:]->temp(event,:,:,:)
level = fhead1[0]->level
leveltmp =  conform_dims(dimsizes(tpl),level,1)
sclhgt = 8000.
R = 290

dthetadz = center_finite_diff_n(tpl*(1000./leveltmp)^0.286,-sclhgt*log(level/1000.),False,0,1)
copy_VarCoords(tpl,dthetadz)
NN = (R*(leveltmp/1000.)^0.286)/sclhgt * dthetadz
copy_VarCoords(tpl,dthetadz)
copy_VarCoords(tpl,NN)

system("rm -f /work1/hyey/data/ecmwf/streamf/dthetadz_NN_cold.b1.nc")
a = addfile("/work1/hyey/data/ecmwf/streamf/dthetadz_NN_cold.b1.nc","c")
a->dthetadz = dthetadz
a->NN = NN 

delete([/a,dthetadz,tpl,fhead1,level,leveltmp,NN/])
flist2 = dir +"Upl/Upl."+cd_calendar(tim1,-1)+".day.nc"
fhead2 = addfiles(flist2,"r")
upl    = fhead2[:]->U_velocity(event,:,::-1,:)
lat = fhead2[0]->lat
delete(fhead2)
print("u")

flist3 = dir +"Vpl/Vpl."+cd_calendar(tim1,-1)+".day.nc"
fhead3 = addfiles(flist3,"r")
vpl    = fhead3[:]->V_velocity(event,:,::-1,:)
lat = fhead3[0]->lat
delete(fhead3)
print("v")

sf = uv2sfvpF(upl,vpl)
sf2 = sf(0,:,:,:,:)
;sf2@_FillValue = tpl@_FillValue
pi = atan(1.0)*4.
d2r = pi/180
a = 6371220.
y = lat*d2r*a
uu = -1*center_finite_diff_n(sf2,y,False,0,2)
copy_VarCoords(upl,sf2)
copy_VarCoords(upl,uu)
;delete([/upl,vpl/])
printVarSummary(sf2)
uu!2 = "lat"
uu&lat = lat
vrt = uv2vrF_Wrap(upl,vpl)  
psi = ilapsF_Wrap(vrt,0)

    system("rm -f /work1/hyey/data/ecmwf/streamf/streamf_cold.b1.nc")
n = addfile("/work1/hyey/data/ecmwf/streamf/streamf_cold.b1.nc","c")
n->sf = sf2(:,:,::-1,:)
n->uu = uu(:,:,::-1,:)
n->ou = upl(:,:,::-1,:)
n->psi = psi(:,:,::-1,:)


