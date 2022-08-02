; Use: Calculate the R.H.S. terms in Quasi-geostrophic omega equation of eq. (6.53) in Holton 4ed using 2nd-order centered difference.
;      1) Q-vector divergence (T42), 2) beta forcing (T42), 3) diabatic forcing (T128).
; : Chanil Park @ Seoul National University
; Last Revision:
load "./ncl_myfunction.ncl"

; Directories of input & output data
dir_in1 = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/levels/Z/value_pr/"
dir_in2 = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/levels/T/value_pr/"
dir_in3 = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/levels/Tphy/value_pr/"
dir_out = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/OmegaEq/RHS_forcing_pr/"


ago = (/"36ago","18ago","00"/)
trunc=True

; Set parameters
rad = 4.0*atan(1.0)/180.0 ; radian factor
a=6.37122e6     ; raidus of earth
omega=7.2921e-5 ; angular veloity of earth
g=9.8 ; gravitational accerelation
R=287 ; universal gas constant
Cp=1004.5 ; specific heat for dry gas under constant pressure
reflat=38


do ago_n = 0,dimsizes(ago)-1
print("-------------- "+ago(ago_n))

;================================================================
; [1] Q-vector forcing (first term)
;================================================================
print("Caculation of Q-vector forcing term")
f_in1 = addfile(dir_in1+"unified."+ago(ago_n)+".Z.nc", "r")
f_in2 = addfile(dir_in2+"unified."+ago(ago_n)+".T.nc", "r")
Gp = f_in1->Z(:,{1000*100:100*100},::-1,:)
T = f_in2->T(:,{1000*100:100*100},::-1,:)
 lon = tofloat(T&lon)
 lat = tofloat(T&lat)
 plev = tofloat(T&lev) ; [Pa]
 time = tofloat(T&time)

if (trunc.eq.True) then
 Tr=42
 spec_harmonics = shaeC(Gp) ; to spherical harmonics domain
 spec_harmonics = tri_trunC(spec_harmonics,Tr) ; truncate
 Gp = shseC(spec_harmonics,dimsizes(lon)) ; return to physical domain
 delete(spec_harmonics)
 spec_harmonics = shaeC(T) ; to spherical harmonics domain
 spec_harmonics = tri_trunC(spec_harmonics,Tr) ; truncate
 T = shseC(spec_harmonics,dimsizes(lon)) ; return to physical domain
 delete(spec_harmonics)
end if

; Q vector divergence
f = 2*omega*sin(lat*rad) ; Coriolis parameter
f0 = 2*omega*sin(reflat*rad) ; Coriolis parameter for reference latitude

result = qvector1(Gp,T,f0,plev,1,0) ; Q-vector using geostrophic wind
delete([/Gp,T/])

Qdiv = result[2]
Qdiv_1 = result[5]
Qdiv_2 = result[8]
Qdiv = where(ismissing(Qdiv),0,Qdiv)
Qdiv_1 = where(ismissing(Qdiv_1),0,Qdiv_1)
Qdiv_2 = where(ismissing(Qdiv_2),0,Qdiv_2)
Qdiv@_FillValue = -999.
Qdiv_1@_FillValue = -999.
Qdiv_2@_FillValue = -999.

Qdiv(:,0,:,:) = 0
Qdiv_1(:,0,:,:) = 0
Qdiv_2(:,0,:,:) = 0

Qdiv = where(abs(Qdiv).gt.1e-5,0,Qdiv)
Qdiv_1 = where(abs(Qdiv_1).gt.1e-5,0,Qdiv_1)
Qdiv_2 = where(abs(Qdiv_2).gt.1e-5,0,Qdiv_2)
out = dir_out+"unified."+ago(ago_n)+".Qforcing.nc"
system("rm -f "+out)
f_out = addfile(out, "c")
f_out->Qfor = Qdiv
f_out->Qparfor = Qdiv_1
f_out->Qperfor = Qdiv_2
delete([/f_in1,f_in2, out, f_out/])
delete([/lon,lat,plev,time/])
delete([/f,f0, result, Qdiv,Qdiv_1,Qdiv_2/])
print("			done") 


;================================================================
; [2] planetary vorticity advection by thermal wind (second term)
;================================================================
print("Caculation of beta effect term forcing")
f_in1 = addfile(dir_in1+"unified."+ago(ago_n)+".Z.nc", "r")
Gp = f_in1->Z(:,{1000*100:100*100},::-1,:)
 lon = tofloat(Gp&lon)
 lat = tofloat(Gp&lat)
 plev = tofloat(Gp&lev) ; [Pa]
 time = tofloat(Gp&time)

if (trunc.eq.True) then
 Tr=42
 spec_harmonics = shaeC(Gp) ; to spherical harmonics domain
 spec_harmonics = tri_trunC(spec_harmonics,Tr) ; truncate
 Gp = shseC(spec_harmonics,dimsizes(lon)) ; return to physical domain
 delete(spec_harmonics)
end if

; beta (df/dy)
f = 2*omega*sin(lat*rad) ; Coriolis parameter
f0 = 2*omega*sin(reflat*rad) ; Coriolis parameter for reference latitude
dy = a*abs(lat(1)-lat(0))*rad
beta = center_finite_diff(f,dy,False,0)
beta0 = beta(minind(lat-reflat))

; meridional geostrophic wind
vg = new(dimsizes(Gp),"float")
dx = new((/dimsizes(lat)/),"float")
do j=0,dimsizes(lat)-1
 dx(j) = a*abs(cos(lat(j)*rad)) * abs(lon(1)-lon(0))*rad
end do
do j=0,dimsizes(lat)-1
  vg(:,:,j:j,:) = 1/f0 * center_finite_diff_n(Gp(:,:,j:j,:),dx(j),True,0,3)
end do

; vertical derivative of the meridional geostrophic wind
dvgdp = center_finite_diff_n(vg,plev,False,2,1)
BETAeffect = f0 * beta0 * dvgdp
copy_VarCoords(Gp,BETAeffect)
BETAeffect = where(ismissing(BETAeffect),0,BETAeffect)
BETAeffect@_FillValue = -999.

BETAeffect = where(abs(BETAeffect).gt.1e-5,0,BETAeffect)
out = dir_out+"unified."+ago(ago_n)+".BETAforcing.nc"
system("rm -f "+out)
f_out = addfile(out, "c")
f_out->betafor = BETAeffect
delete([/f_in1, out, f_out/])
delete([/Gp,lon,lat,plev,time/])
delete([/f,f0, dy,beta, vg,dx, dvgdp, BETAeffect/])
print("			done")
 
;================================================================
; [3] diabatic forcing (third temr)
;================================================================
print("Caculation of diabatic forcing term")
f_in1 = addfile(dir_in3+"unified."+ago(ago_n)+".Tphy.nc", "r") ; Tphy data is 6-hour forward averate (ex. 00: 00-06 average)
Tphy = f_in1->Tphy(:,{1000*100:100*100},::-1,:)
 lon = tofloat(Tphy&lon)
 lat = tofloat(Tphy&lat)
 plev = tofloat(Tphy&lev) ; [Pa]
 time = tofloat(Tphy&time)

if (trunc.eq.True) then
 Tr=128
 Tphy = where(ismissing(Tphy),0,Tphy)
 spec_harmonics = shaeC(Tphy) ; to spherical harmonics domain
 spec_harmonics = tri_trunC(spec_harmonics,Tr) ; truncate
 Tphy = shseC(spec_harmonics,dimsizes(lon)) ; return to physical domain
 delete(spec_harmonics)
end if

Tphy = Tphy/(60*60*6) ; K/6hours -> K/s

; Laplacian of total physics tendency (K/s)
dx = new((/dimsizes(lat)/),"float")
do j=0,dimsizes(lat)-1
 dx(j) = a*abs(cos(lat(j)*rad)) * abs(lon(1)-lon(0))*rad
end do
dy = a*abs(lat(1)-lat(0))*rad

d2Tphydx2 = Tphy*0
do j=0,dimsizes(lat)-1
 d2Tphydx2(:,:,j:j,:) = center_finite_diff_n(Tphy(:,:,j:j,:),dx(j),False,2,3)
end do
do j=0,dimsizes(lat)-1
 d2Tphydx2(:,:,j:j,:) = center_finite_diff_n(d2Tphydx2(:,:,j:j,:),dx(j),False,2,3)
end do

d2Tphydy2 = Tphy*0
d2Tphydy2 = center_finite_diff_n(Tphy,dy,False,2,2)
d2Tphydy2 = center_finite_diff_n(d2Tphydy2,dy,False,2,2)

LapTphy = Tphy*0
LapTphy = (d2Tphydx2 + d2Tphydy2)

; final form of diabatic forcing term
do k=0,dimsizes(plev)-1
 LapTphy(:,k,:,:) = (-1/plev(k)) * (R/Cp) * LapTphy(:,k,:,:)*Cp
end do
copy_VarCoords(Tphy,LapTphy)
LapTphy = where(ismissing(LapTphy),0,LapTphy)
LapTphy@_FillValue = -999.

LapTphy = where(abs(LapTphy).gt.1e-5,0,LapTphy)
out = dir_out+"unified."+ago(ago_n)+".DIAforcing.nc"
system("rm -f "+out)
f_out = addfile(out, "c")
f_out->diafor = LapTphy
delete([/f_in1, out, f_out/])
delete([/Tphy,lon,lat,plev,time/])
delete([/dx,dy, d2Tphydx2,d2Tphydy2, LapTphy/])
print("                         done")


end do
