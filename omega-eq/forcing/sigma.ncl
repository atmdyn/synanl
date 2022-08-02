; Directories of input & output data
dir_in1 = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/levels/T/value_pr/"
dir_out = "/home/storage-eddy/cipark/Processed/Myresearch/Heavyrain2/OmegaEq/RHS_forcing_pr/"

ago = (/"36ago","18ago","00"/)

; Set parameters
rad = 4.0*atan(1.0)/180.0 ; radian factor
a=6.37122e6     ; raidus of earth
omega=7.2921e-5 ; angular veloity of earth
g=9.8 ; gravitational accerelation
R=287 ; dry gas constant
Cp=1004.5 ; specific heat for dry gas under constant pressure

do ago_n = 0,dimsizes(ago)-1

f_in1 = addfile(dir_in1+"unified"+"."+ago(ago_n)+".T.nc", "r")
T = f_in1->T(:,{1000*100:50*100},{10:65},{85:170})
 lon = tofloat(T&lon)
 lat = tofloat(T&lat) 
 plev = tofloat(T&lev) ; [Pa]
 time = tofloat(T&time)

; areal mean temperature for air density and potential temperature
T2 = T(:,:,{28:48},{110:145})

lat2 = tofloat(T2&lat)
T_area = wgt_areaave(T2, cos(lat2*rad), 1.0, 1)

; air density and potential temperature
rho = T_area*0 ; density of air
do k=0,dimsizes(plev)-1
 rho(:,k) = plev(k) / ( T_area(:,k)*R )
end do
PT = pot_temp(plev,T_area,1,False)


; static stability
dPTdp = center_finite_diff_n(PT,plev,False,2,1)
sigma = dPTdp*0
do t=0,dimsizes(time)-1
do k=0,dimsizes(plev)-1
 sigma(t,k) = -(1/rho(t,k)) * (1/PT(t,k)) * dPTdp(t,k)
end do
end do
sigma2 = T*0
sigma2(:,:,0,0) = sigma
copy_VarCoords(T,sigma2)

sigma2 = where(ismissing(sigma2),0,sigma2)
sigma2@_FillValue = -999.

out = dir_out+"unified"+"."+ago(ago_n)+".sigma.nc"
system("rm -f "+out)
f_out = addfile(out, "c")
f_out->sigma = sigma2(:,:,0:0,0:0)

delete([/f_in1, T, lon,lat,plev,time, T2, lat2, T_area, rho,PT, dPTdp,sigma,sigma2, out,f_out/])
end do

print("                         done")
