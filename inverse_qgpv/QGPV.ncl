fi = addfile("/work1/hyey/data/ecmwf/streamf/qgpv/v2/streamf_cold.b1.nc","r")
sf = fi->sf
day = fi->day
lat = fi->lat
lon = fi->lon
lev = fi->level
nday = dimsizes(day)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
nlev = dimsizes(lev)
lev4 = conform_dims(dimsizes(sf),lev,1)

fi2 = addfile("/work1/hyey/data/ecmwf/streamf/dthetadz_NN_cold.b1.nc","r")
NN = fi2->NN
NN := dim_avg_n(dim_avg_n(NN,2),2)
NN := conform_dims(dimsizes(sf),NN(0,:),1)
NN = where(NN .gt. 0, NN, 1e-5)
R = 287. ; J/kg*K
H = 7500.

f = coriolis_param(lat)
f := conform_dims(dimsizes(sf),f,2)

f0 = 10^-4    
pi = atan(1.0)*4.
d2r = pi/180
a = 6371220.
y = lat*d2r*a
print(y)
    exit
A = f^2/NN
sf_x = new((/nday,nlev,nlat,nlon/),float)
sf_xx = sf_x
do i = 0, dimsizes(lat)-1
  x = lon*d2r*a*cos(lat(i)*d2r)

  sf_x(:,:,i,:) = center_finite_diff_n(sf(:,:,i,:),x,True,0,2)    
  sf_xx(:,:,i,:) = center_finite_diff_n(sf_x(:,:,i,:),x,True,0,2)
end do

sf_yy = sf_xx
sf_yy = 0 
j = ispan(1,dimsizes(lat)-2,1)
dy = y(2)-y(1)
sf_yy(:,:,j,:) = (sf(:,:,j+1,:)-2*sf(:,:,j,:)+sf(:,:,j-1,:))/(dy^2)


z = -H*log(lev/1000.)
sf_z = center_finite_diff_n(sf,z,False,0,1)
sf_zz = center_finite_diff_n(lev4*(f0^2/NN*sf_z),z,False,0,1)/lev4

printVarSummary(sf_xx)
printVarSummary(sf_yy)
printVarSummary(sf_zz)
qgpv =  sf_xx + sf_yy + sf_zz




copy_VarCoords(sf,qgpv)
system("rm -f test.nc")
nfi = addfile("test.nc","c")
nfi->qgpv = qgpv
nfi->mqgpv = dim_avg_n_Wrap(qgpv,0)
;system("rm -f /work1/hyey/data/ecmwf/streamf/qgpv/v2/qgpv_cold.b1.nc")
;nfi = addfile("/work1/hyey/data/ecmwf/streamf/qgpv/v2/qgpv_cold.b1.nc","c")
;nfi->qgpv = qgpv
;nfi->mqgpv = dim_avg_n_Wrap(qgpv,0)
;nfi->sf_xx =sf_xx
;nfi->sf_yy =sf_yy
;nfi->sf_zz =sf_zz
exit
kc1 = asciiread("/home/hyey/coldsurge/2019/synop/first_cold/clust/runavg/new_cl/kc1.txt",28,"integer")
kc2 = asciiread("/home/hyey/coldsurge/2019/synop/first_cold/clust/runavg/new_cl/kc2.txt",33,"integer")

qgpv1 = dim_avg_n_Wrap(qgpv(kc1,:,:,:),0)    
qgpv2 = dim_avg_n_Wrap(qgpv(kc2,:,:,:),0)    
    printVarSummary(qgpv1)

system("rm -f /work1/hyey/data/ecmwf/streamf/qgpv_kc1.b1.nc")
b = addfile("/work1/hyey/data/ecmwf/streamf/qgpv_kc1.b1.nc","c")
b->qgpv = qgpv(kc1,:,:,:)

system("rm -f /work1/hyey/data/ecmwf/streamf/mqgpv_kc1.b1.nc")
b1 = addfile("/work1/hyey/data/ecmwf/streamf/mqgpv_kc1.b1.nc","c")
b1->qgpv = qgpv1

system("rm -f /work1/hyey/data/ecmwf/streamf/qgpv_kc2.b1.nc")
c = addfile("/work1/hyey/data/ecmwf/streamf/qgpv_kc2.b1.nc","c")
c->qgpv = qgpv(kc2,:,:,:)

system("rm -f /work1/hyey/data/ecmwf/streamf/mqgpv_kc2.b1.nc")
c1 = addfile("/work1/hyey/data/ecmwf/streamf/mqgpv_kc2.b1.nc","c")
c1->qgpv = qgpv2
