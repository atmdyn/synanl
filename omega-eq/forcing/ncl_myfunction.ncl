;		[CUSTOMIZED NCL FUNCTIONS FOR DYNAMIC-SYNOPTIC METEOROLOGY ANALYSIS]
; 				  Last revision: 2021/04/11
;
; Function list: fronto1, fronto2, qvector1, vor_CurSh, TNWAF_2d
;
; Two types of horizontal differentiation are supported:
;   1) finite central differnce (grad_latlon_cfd)
;
;   2) difference using spherical harmonics (gradsf).
;       - When using spherical harmonics (using 'gradsf'), variables should have globaly-distributed both in longitude and latitude, and the latitude should be in S->N order.
;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: fronto1
; Use: Compute 2-D frontogenesis function by Petterssen (1936,1956)
; Output: frontogenesis Fn by 1) deformation 2) divergence 
;			Fs by 1) vorticity   2) deformation
; Developer: Chanil Park @ Seoul National University
; Last revision: 2019.11.25
;
; Arguments:
; theta : (equivalent) potential temperature (K, 4-d or 3-d data)
; u     : zonal component of flow (m/s, 4-d or 3-d data)
; v     : meridional component of flow (m/s, 4-d or 3-d data)
; opt   : 0 for finite central difference, 1 for use of spherical harmonics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("fronto1")
function fronto1 (theta:float, u:float, v:float, opt:integer)
local dvdx, dvdy, dudy, dudx, dthdy, dthdx, mag_grad_theta, vor, div, Est, Esh, def, alpha, delta, beta
begin

 if (opt.eq.0) then
  if ((dimsizes(dimsizes(v))).eq.4) then
   u!2="lat"
   u!3="lon"
  else if ((dimsizes(dimsizes(v))).eq.3) then
   u!1="lat"
   u!2="lon"
  end if
  end if
  lat=tofloat(u&lat)
  lon=tofloat(u&lon)
  asdf = grad_latlon_cfd(u,lat,lon,False,False)
  dudx = asdf[1]
  dudy = asdf[0]
  asdf = grad_latlon_cfd(v,lat,lon,False,False)
  dvdx = asdf[1]
  dvdy = asdf[0]
  asdf = grad_latlon_cfd(theta,lat,lon,False,False)
  dthdx = asdf[1]
  dthdy = asdf[0]
 else if (opt.eq.1) then
  dvdx = v
  dvdy = v
  dudy = u
  dudx = u
  dthdx = theta
  dthdy = theta
  gradsf(u,dudx,dudy)
  gradsf(v,dvdx,dvdy)
  gradsf(theta,dthdx,dthdy)
 end if
 end if


mag_grad_theta = (dthdx^2.+dthdy^2.)^0.5	; temperature gradient magnitude
vor = dvdx-dudy					; (relative) vorticity
div = dudx+dvdy 				; divergence
Est = dudx-dvdy 				; stretching deformation
Esh = dvdx+dudy 				; shearing deformation
Est@_FillValue = -999.
Esh@_FillValue = -999.
def = (Esh^2.+Est^2.)^0.5			; total deformation
Est = where(Est.eq.0,Est@_FillValue,Est)
delta = atan(Esh/Est)/2.

dthdy@_FillValue = -999.
dthdy = where(dthdy.eq.0,dthdy@_FillValue,dthdy)
alpha = atan(-1.*dthdx/dthdy)
beta = delta-alpha

; Terms related to magnitude of temperature gradient
Fn_def = 0.5*mag_grad_theta*def*cos(2.*beta)*86400.*100000. ; unit: (K/100km)/day
Fn_div = -0.5*mag_grad_theta*div*86400.*100000.

Fn_def@_FillValue = -999.
Fn_def = where(ismissing(Est),Fn_def@_FillValue,Fn_def)
copy_VarCoords(u,Fn_def)
Fn_def@units = "(K/100km)/day"
Fn_def@long_name = "Petterssen frontogenesis (Fn) by deformation"
Fn_div@_FillValue = -999.
Fn_div = where(ismissing(Est),Fn_div@_FillValue,Fn_div)
copy_VarCoords(u,Fn_div)
Fn_div@units = "(K/100km)/day"
Fn_div@long_name = "Petterssen frontogenesis (Fn) by divergence"


; Terms related to direction of temperature gradient
Fs_def = 0.5*mag_grad_theta*def*sin(2.*beta)*86400.*100000. ; unit: (K/100km)/day
Fs_vor = 0.5*mag_grad_theta*vor*86400.*100000.

Fs_def@_FillValue = -999.
Fs_def = where(ismissing(Est),Fs_def@_FillValue,Fs_def)
copy_VarCoords(u,Fs_def)
Fs_def@units = "(K/100km)/day"
Fs_def@long_name = "Petterssen frontogenesis (Fs) by deformation"

Fs_vor@_FillValue = -999.
Fs_vor = where(ismissing(Est),Fs_vor@_FillValue,Fs_vor)
copy_VarCoords(u,Fs_vor)
Fs_vor@units = "(K/100km)/day"
Fs_vor@long_name = "Petterssen frontogenesis (Fs) by vorticity"


return([/Fn_def,Fn_div, Fs_def,Fs_vor/])
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: front2
; Use: Compute full 3-D frontogenesis function by Ninomiya (1984; kinematic function)
; Output: frontogenesis(or lysis) by
;         1) diabatic process 2) divergence 3) deformation 4) tilting
; Caution: Latitudinal grid should be in South-to-North order (advect_variable)
; Developer: Chanil Park @ Seoul National University
; Last revision: 2019.03.16
;
; Arguments:
; theta : potential temperature (K, 4-d or 3-d data)
; u  : zonal component of flow (m/s, 4-d or 3-d data)
; v  : meridional component of flow (m/s, 4-d or 3-d data)
; w  : vertical velocity (Pa/s, 4-d or 3-d data)
; dt  : delta time (s, used for time differentiation)
; plev  : 1-D vertical p-level array (Pa)
; opt   : 0 for finite central difference, 1 for use of spherical harmonics
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("fronto2")
function fronto2 (theta:float, u:float, v:float, w:float, dt:numeric, plev:float, opt:integer)
local dvdx,dvdy,dudy,dudx,dwdx,dwdy,dthdy,dthdx, denominator, div,Est,Esh, th_adv,dthdt,DthDt, dDthDtdx,dDthDtdy,dthdp, dp, FG1,FG2,FG3,FG4
begin
; Preparation
 if (opt.eq.0) then
  if ((dimsizes(dimsizes(v))).eq.4) then
   u!2="lat"
   u!3="lon"
  else if ((dimsizes(dimsizes(v))).eq.3) then
   u!1="lat"
   u!2="lon"
  end if
  end if
  lat=tofloat(u&lat)
  lon=tofloat(u&lon)
  asdf = grad_latlon_cfd(u,lat,lon,False,False)
  dudx = asdf[1]
  dudy = asdf[0]
  asdf = grad_latlon_cfd(v,lat,lon,False,False)
  dvdx = asdf[1]
  dvdy = asdf[0]
  asdf = grad_latlon_cfd(w,lat,lon,False,False)
  dwdx = asdf[1]
  dwdy = asdf[0]
  asdf = grad_latlon_cfd(theta,lat,lon,False,False)
  dthdx = asdf[1]
  dthdy = asdf[0]
 else if (opt.eq.1) then
  dvdx = v
  dvdy = v
  dudy = u
  dudx = u
  dwdy = w
  dwdx = w
  dthdx = theta
  dthdy = theta
  gradsf(u,dudx,dudy)
  gradsf(v,dvdx,dvdy)
  gradsf(w,dwdx,dwdy)
  gradsf(theta,dthdx,dthdy)
 end if
 end if

asdf2 = (dthdx^2.+dthdy^2.)^0.5
asdf2@_FillValue = -999.
asdf2 = where(asdf2.eq.0,asdf2@_FillValue,asdf2)
denominator = 1.0 / asdf2
div = dudx+dvdy         ; divergence
Est = dudx-dvdy         ; stretching deformation
Esh = dvdx+dudy         ; shearing deformation

th_adv = advect_variable(u,v,theta,0,"not important","not important",0)
dthdt = center_finite_diff_n(theta,dt,False,2,0)
DthDt = th_adv + dthdt
copy_VarCoords(u,DthDt)

if (opt.eq.0) then
 asdf = grad_latlon_cfd(DthDt,lat,lon,False,False)
 dDthDtdx = asdf[1]
 dDthDtdy = asdf[0]
else if(opt.eq.1) then
 dDthDtdx = theta
 dDthDtdy = theta
 gradsf(DthDt,dDthDtdx,dDthDtdy)
end if
end if

dthdp = center_finite_diff_n(theta,plev,False,2,1)
;np=dimsizes(plev)
;dp = new(np,"float")
;dp(0)      = exp((log(plev(0))+log(plev(1)))*0.5)              ; set bottom
;dp(np-1)   = exp((log(plev(np-1))+log(plev(np-2)))*0.5)        ; set top    
;dp(1:np-2) = exp((log(plev(0:np-3)) + log(plev(2:np-1)))*0.5)  ; mid points
;dthdp = center_finite_diff_n(theta,dp,False,2,1)


; Frontogenesis terms
FG1 = denominator * ( dthdx*dDthDtdx + dthdy*dDthDtdy )
FG2 = -denominator/2 * ( dthdx^2.+dthdy^2. ) * div
FG3 = -denominator * ( (dthdx^2.-dthdy^2.)/2.*Est ) + dthdx*dthdy*Esh
FG4 = -denominator * dthdp*(dwdx*dthdx+dwdy*dthdy)
FG1 = FG1*86400.*100000.
FG2 = FG2*86400.*100000. ; unit: (K/100km)/day
FG3 = FG3*86400.*100000.
FG4 = FG4*86400.*100000.


copy_VarCoords(u,FG1)
FG1@_FillValue = -999.
FG1@units = "(K/100km)/day"
FG1@long_name = "frontogenesis by diabatic processes"
copy_VarCoords(u,FG2)
FG2@_FillValue = -999.
FG2@units = "(K/100km)/day"
FG2@long_name = "frontogenesis by divergence"
copy_VarCoords(u,FG3)
FG3@_FillValue = -999.
FG3@units = "(K/100km)/day"
FG3@long_name = "frontogenesis by deformation"
copy_VarCoords(u,FG4)
FG4@_FillValue = -999.
FG4@units = "(K/100km)/day"
FG4@long_name = "frontogenesis by tilting"

return([/FG1,FG2,FG3,FG4/])
end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: qvector1
; Use: Calculate Q vector (based on eq. (6.54) in Holton 4ed) and Q vector forcing (Q-vector divergence muliplied by -2)
;      See Eq. (6.54) in Holton 4ed for the Q vector formula in Cartesian coordinate
; Output: Depends on 'opt1' (see below)
; Caution: In the lowest level, all values are undifined.
;         (for example, if you use 1000-250 hPa levels data, result at 1000 hPa is unfined, but I don't know the reason...)
; Developer: Chanil Park @ Seoul National University
; Last Revision: 2019.08.24
;
; Arguments:
; Gp  : geopotential (m2/s2, 4-d data)
; T  : temperature (K, 4-d data)
; f  : Coriolis parameter for each latitude (same dimension with latitude dimension) or at reference latitude (scalar) for QG approximation
; plev  : 1-D vertical p-level array (Pa)
; opt1  : 0 or 1; see below
; opt2  : 0 for finite central difference, 1 for use of spherical harmonics
;
; opt1 0: Q vector and its forcing.
;            output: ([/Qx,Qy,Qdiv/])
;      1: additionally, two independent components of Q vector (shearwise;Qs vs transverse;Qn) and their forcings.
;            output: ([/Qx,Qy,Qdiv, Qx_s,Qy_s,Qdiv_s, Qx_n,Qy_n,Qdiv_n/])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("qvector1")
function qvector1 (Gp:float, T:float, f:float, plev:float, opt1:integer, opt2:integer)
local rad,a,omega,g,R, dGpdx,dGpdy,ug,vg,dTdx,dTdy, Qx,Qy,Qdiv, denominator, denominator2, Qx_1,Qy_1,dQxdx_1,dQxdy_1,dQiv_1, Qx_2,Qy_2,dQxdx_2,dQxdy_2,dQiv_2
begin
 R=287 ; dry gas constant

; Calculate geostrophic wind 
 if (opt2.eq.0) then
  Gp!2="lat"
  Gp!3="lon"
  lat=tofloat(Gp&lat)
  lon=tofloat(Gp&lon)
  asdf = grad_latlon_cfd(Gp,lat,lon,False,False)
  dGpdx = asdf[1]
  dGpdy = asdf[0] 
 else if (opt2.eq.1) then
  dGpdx = Gp
  dGpdy = Gp
  gradsf(Gp,dGpdx,dGpdy)
 end if  
 end if

 fdims = dimsizes(f)
 f = where(f.eq.0,1,f) ; coriolis parameter at equator 
 if (fdims.eq.1) then	; if using reference latitude for Coriolis parameter
	 ug = -1/f * dGpdy
	 vg = 1/f * dGpdx
 else if(fdims.ge.2) then ; if using latitude-dependent Coriolis parameter
	ug = new(dimsizes(Gp),"float")
	vg = new(dimsizes(Gp),"float")
	do j=0,fdims-1
	 ug(:,:,j,:) = -1/f(j) * dGpdy(:,:,j,:)
     vg(:,:,j,:) = 1/f(j) * dGpdx(:,:,j,:)
	end do
 end if
 end if	
 copy_VarCoords(Gp,ug)
 copy_VarCoords(Gp,vg)


; Calculate Q-vector
 if (opt2.eq.0) then
  asdf = grad_latlon_cfd(ug,lat,lon,False,False)
  dugdx = asdf[1]
  dugdy = asdf[0]
  asdf = grad_latlon_cfd(vg,lat,lon,False,False)
  dvgdx = asdf[1]
  dvgdy = asdf[0]
  asdf = grad_latlon_cfd(T,lat,lon,False,False)
  dTdx = asdf[1]
  dTdy = asdf[0]
 else if (opt2.eq.1) then
  dugdx=Gp
  dugdy=Gp
  dvgdx=Gp
  dvgdy=Gp
  dTdx=Gp
  dTdy=Gp
  gradsf(ug,dugdx,dugdy)
  gradsf(vg,dvgdx,dvgdy)
  gradsf(T,dTdx,dTdy)
 end if
 end if

 Qx = new(dimsizes(Gp),"float")
 Qy = new(dimsizes(Gp),"float") 
 do k=1,dimsizes(plev)-1
  Qx(:,k,:,:) = -1* R/plev(k) * (dugdx(:,k,:,:)*dTdx(:,k,:,:) + dvgdx(:,k,:,:)*dTdy(:,k,:,:))
  Qy(:,k,:,:) = -1* R/plev(k) * (dugdy(:,k,:,:)*dTdx(:,k,:,:) + dvgdy(:,k,:,:)*dTdy(:,k,:,:))
 end do
 copy_VarCoords(T,Qx)
 copy_VarCoords(T,Qy)

 if (opt2.eq.0) then
  Qdiv = uv2dv_cfd(Qx,Qy,lat,lon,2)
 else if (opt2.eq.1) then
  Qdiv = Gp
  uv2dvf(Qx,Qy,Qdiv)
 end if
 end if
 Qdiv = -2*Qdiv
 copy_VarCoords(T,Qdiv)

 Qx@_FillValue = -999.
 Qy@_FillValue = -999.
 Qdiv@_FillValue = -999.
 Qx@long_name = "Q vector x component"
 Qy@long_name = "Q vector y component"
 Qdiv@long_name = "Q vector divergence"


 if (opt1.eq.1) then
   denominator = dTdx^2 + dTdy^2
   denominator@_FillValue = -999.
   denominator = where(denominator.eq.0, denominator@_FillValue, denominator) ; just for case of being devided by 0

   ; Trensverse component
   Qx_1 = ( Qx*dTdx^2 + Qy*dTdx*dTdy ) / denominator
   Qy_1 = ( Qx*dTdx*dTdy + Qy*dTdy^2 ) / denominator
   copy_VarCoords(T,Qx_1)
   copy_VarCoords(T,Qy_1)
   if (opt2.eq.0) then
    Qdiv_1 = uv2dv_cfd(Qx_1,Qy_1,lat,lon,2)
   else if (opt2.eq.1) then
    Qdiv_1 = Gp
    uv2dvf(Qx_1,Qy_1,Qdiv_1)
   end if
   end if
   Qdiv_1 = -2*Qdiv_1
   copy_VarCoords(T,Qdiv_1)
   Qx_1@_FillValue = -999.
   Qy_1@_FillValue = -999.
   Qdiv_1@_FillValue = -999.
   Qx_1@long_name = "Q vector x component (frontogenetic, transverse)"
   Qy_1@long_name = "Q vector y component (frontogenetic, transverse)"
   Qdiv_1@long_name = "Q vector divergence (frontogenetic, transverse)"
   ; Shearwise component
   Qx_2 = ( -Qy*dTdx*dTdy + Qx*dTdy^2 ) / denominator
   Qy_2 = ( Qy*dTdx^2 - Qx*dTdx*dTdy ) / denominator
   copy_VarCoords(T,Qx_2)
   copy_VarCoords(T,Qy_2)
   if (opt2.eq.0) then
    Qdiv_2 = uv2dv_cfd(Qx_2,Qy_2,lat,lon,2)
   else if (opt2.eq.1) then
    Qdiv_2 = Gp
    uv2dvf(Qx_2,Qy_2,Qdiv_2)
   end if
   end if
   Qdiv_2 = -2*Qdiv_2
   copy_VarCoords(T,Qdiv_2)
   Qx_2@_FillValue = -999.
   Qy_2@_FillValue = -999.
   Qdiv_2@_FillValue = -999.
   Qx_2@long_name = "Q vector x component (rotating, shearwise)"
   Qy_2@long_name = "Q vector y component (rotating, shearwise)"
   Qdiv_2@long_name = "Q vector divergence (rotating, shearwise)"

   return ([/Qx,Qy,Qdiv, Qx_1,Qy_1,Qdiv_1, Qx_2,Qy_2,Qdiv_2/])

 else if (opt1.eq.0) then
   return ([/Qx,Qy,Qdiv/])
 end if
 end if

end


;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: qvector2
; Use: Calculate shearwise (Qs) and transverse Q vectors (Qn) and their forcings in kinematic sense.
;      The shearwise Q-vector comes from geostrophic relative vorticity (VR) and geostrohpic deformation (DR),
;      and the transverse Q-vector comes from geostrophic deformation (Keyser et al. 1988; Martin 1999).
;      [note: Each component of Qs forcings is identical to VR or DR advection by thermal wind, see Martin (1999)]
; Output: Depends on 'opt1' (see below)
; Caution:
; Developer: Chanil Park @ Seoul National University
; Last Revision:
;
; Arguments:
; Gp  : geopotential (m2/s2, 4-d data)
; theta  : potential temperature (K, 4-d data)
; f  : Coriolis parameter for each latitude (same dimension with latitude dimension) or at reference latitude (scalar) for QG approximation
; plev  : 1-D vertical p-level array (Pa)
; opt1  : 0 or 1 or 2; see below
; opt2  : 0 for finite central difference, 1 for use of spherical harmonics
;
; opt1 0: x and y components of Qs vector by geostrophic relative vorticity and its forcing
;		output: ([/Qx_s_VR, Qy_s_VR, Q_s_VR_div/])
;      1: x and y components of Qs vector by geostrophic deformation and its forcing
;               output: ([/Qx_s_DR, Qy_s_DR, Q_s_DR_div/])
;      2: x and y components of Qn vector by geostrophic deformation and its forcing
;               output: ([/Qx_n_DR, Qy_n_DR, Q_n_DR_div/])
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;

;
; under construction ...
;



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: vor_CurSh
; Use: Compute vorticity by curvature and vorticity by shear (eq.1 and eq.2 in 'Lackmann et al. (1997)' or see 'Bell and Keyser (1993)')
; Output: Vorticity by
;         1) curvature 2) shear
; Developer: Chanil Park @ Seoul National University
; Last revision: 2019.04.18
;
; Arguments:
; u  : zonal component of flow (m/s, 4-d or 3-d data)
; v  : meridional component of flow (m/s, 4-d or 3-d data)
; lat  : one dimensional latitude array (south-to-north)
; opt1 : 0 for finite central difference, 1 for use of spherical harmonics
; opt2 : 0 for including earth curvature and planetary vorticity in curvature vorticity, 1 for not including 
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("vor_CurSh")
function vor_CurSh (u:float, v:float, lat:float, opt1:integer, opt2:integer)
local rad,a,omega, f,Cur_earth, wind, dvdx, dvdy, dudy, dudx, vor_cur, vor_sh
begin

rad = 4.0*atan(1.0)/180.0 ; radian factor
a=6.37122e6     ; raidus of earth
omega=7.2921e-5 ; angular veloity of earth

 if (opt1.eq.0) then
  if ((dimsizes(dimsizes(v))).eq.4) then
   u!2="lat"
   u!3="lon"
  else if ((dimsizes(dimsizes(v))).eq.3) then
   u!1="lat"
   u!2="lon"
  end if
  end if
  lat=tofloat(u&lat)
  lon=tofloat(u&lon)
  asdf = grad_latlon_cfd(u,lat,lon,False,False)
  dudx = asdf[1]
  dudy = asdf[0]
  asdf = grad_latlon_cfd(v,lat,lon,False,False)
  dvdx = asdf[1]
  dvdy = asdf[0]
 else if (opt1.eq.1) then
  dvdx = v
  dvdy = v
  dudy = u
  dudx = u
  gradsf(u,dudx,dudy)
  gradsf(v,dvdx,dvdy)
 end if
 end if

wind = u^2+v^2
wind = where(wind.eq.0,1e18,wind) ; for case of wind=0 (cannot divide!)
; Shear vorticity
vor_sh = 1/wind * (u*v*(dudx-dvdy) + v^2*dvdx - u^2*dudy)

; Curvature vorticity
vor_cur = 1/wind * (-u*v*(dudx-dvdy) + u^2*dvdx - v^2*dudy )

if(opt2.eq.0) then
 f=v
 f=0
 Cur_earth = v
 Cur_earth = 0
   if (dimsizes(dimsizes(v)).eq.4) then ; time, lev, lat, lon
      do j=0,dimsizes(lat)-1
         f(:,:,j,:) = where(lat(j).eq.0,1,2*omega*sin(lat(j)*rad)) ; Coriolis parameter
         Cur_earth(:,:,j,:) = u(:,:,j,:)*tan(lat(j)*rad)/a ; Curvature by spherical earth
      end do
   else if(dimsizes(dimsizes(v)).eq.3) then ; time, lat, lon
      do j=0,dimsizes(lat)-1
         f(:,j,:) = where(lat(j).eq.0,1,2*omega*sin(lat(j)*rad)) ; Coriolis parameter
         Cur_earth(:,j,:) = u(:,j,:)*tan(lat(j)*rad)/a ; Curvature by spherical earth
      end do
   end if
   end if
 vor_cur = vor_cur + Cur_earth + f
else if(opt2.eq.1) then
end if
end if


copy_VarCoords(v,vor_cur)
vor_cur@_FillValue = -999.
if (opt2.eq.0) then
	vor_cur@long_name = "vorticity by curvature"
else if(opt2.eq.1) then
	vor_cur@long_name = "vorticity by curvature without earth curvature and planetary vorticity"
end if
end if

copy_VarCoords(v,vor_sh)
vor_sh@_FillValue = -999.
vor_sh@long_name = "vorticity by shear"
return ([/vor_cur, vor_sh/])

end



;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
; Function: TNWAF_2d
; Use: Compute 2-D Takaya-Nakamura wave activity flux (eq. 38 in Takaya and Nakamura 2001)
; Output: i and j components of Takaya-Nakamura wave activity flux
; Developer: Chanil Park @ Seoul National University
; Last revision: 2020.11.11
;
; Note: 1) TNWAF_2d is applied to 3-D (nlev x nlat x nlon) or 2-D (nlat x nlon) data.
;
; Arguments:
; psia : stream function anomaly (geostrophic or full)
; U  : background zonal wind (m/s, 3-d or 2-d data, geostrophic or full)
; V  : background meridional wind (m/s, 3-d or 2-d data, geostrophic or full)
; plev  : vertical pressure level (1-D or scalar; Pa)
; opt1 : 0 for finite central difference, 1 for use of spherical harmonics
; opt2 : 0 for omitting pressure factor, 1 for retainment
;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
undef("TNWAF_2d")
function TNWAF_2d (psia:float, U:float, V:float, plev:float, opt1:integer, opt2:integer)
local rad,a, psia,U,V,plev, dpsiadx,dpsiady, d2psiadx2, d2psiadxdy, d2psiady2, WAFx,WAFy

begin

rad = 4.0*atan(1.0)/180.0 ; radian factor
a=6.37122e6     ; raidus of earth


 if (opt1.eq.0) then
  if ((dimsizes(dimsizes(psia))).eq.3) then
   psia!1="lat"
   psia!2="lon"
  else if ((dimsizes(dimsizes(psia))).eq.2) then
   psia!0="lat"
   psia!1="lon"
  end if
  end if
  lat=tofloat(psia&lat)
  lon=tofloat(psia&lon)
  asdf = grad_latlon_cfd(psia,lat,lon,False,False)
  dpsiadx = asdf[1]
  dpsiady = asdf[0]
 else if (opt1.eq.1) then
  dpsiadx = psia
  dpsiady = psia
  gradsf(psia,dpsiadx,dpsiady)
 end if
 end if

 if (opt1.eq.0) then
  asdf = grad_latlon_cfd(dpsiadx,lat,lon,False,False)
  d2psiadx2  = asdf[1]
  d2psiadxdy = asdf[0]
 else if (opt1.eq.1) then
  d2psiadx2 = psia
  d2psiadxdy = psia
  gradsf(dpsiadx,d2psiadx2,d2psiadxdy)
 end if
 end if

 if (opt1.eq.0) then
  asdf = grad_latlon_cfd(dpsiady,lat,lon,False,False)
  d2psiady2 = asdf[0]
 else if (opt1.eq.1) then
  d2psiadxdy = psia
  d2psiady2 = psia
  gradsf(dpsiady,d2psiadxdy,d2psiady2)
 end if
 end if


factor = psia*0
if (opt2.eq.0) then
	do j=0,dimsizes(lat)-1
	 if ((dimsizes(dimsizes(psia))).eq.3) then
	  factor(:,j,:) = cos(lat(j)*rad)/(2*sqrt(U(:,j,:)^2+V(:,j,:)^2))
	 else if ((dimsizes(dimsizes(psia))).eq.2) then
	  factor(j,:) = cos(lat(j)*rad)/(2*sqrt(U(j,:)^2+V(j,:)^2))
	 end if
	 end if
	end do

else if (opt2.eq.1) then
	do j=0,dimsizes(lat)-1
	 if ((dimsizes(dimsizes(psia))).eq.3) then
	  do k=0,dimsizes(plev)-1
	   factor(k,j,:) = plev*cos(lat(j)*rad) / (2*sqrt(U(k,j,:)^2+V(k,j,:)^2))
	  end do

	 else if ((dimsizes(dimsizes(psia))).eq.2) then
	  factor(j,:) = plev*cos(lat(j)*rad) / (2*sqrt(U(j,:)^2+V(j,:)^2))
	 end if
	 end if
	end do
end if
end if

WAFx = factor* ( U*( dpsiadx^2 - psia*d2psiadx2 ) + V*( dpsiadx*dpsiady - psia*d2psiadxdy ) )
WAFy = factor* ( U*( dpsiadx*dpsiady - psia*d2psiadxdy ) + V* ( dpsiady^2 - psia*d2psiady2 ) )

copy_VarCoords(psia,WAFx)
WAFx@_FillValue = -999.
WAFx@long_name = "zonal component of Takaya-Nakamura wave activity flux"
copy_VarCoords(psia,WAFy)
WAFy@_FillValue = -999.
WAFy@long_name = "meridional component of Takaya-Nakamura wave activity flux"

return ([/WAFx, WAFy/])

end
