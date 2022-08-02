
; 2-2.Function : qvector1/*{{{*/
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
;--------------------------------------------------------------------------------------------
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
 if (fdims.eq.1) then   ; if using reference latitude for Coriolis parameter
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

end;/*}}}*/

; 2-3.Function : OmegaEqSolver/*{{{*/
;
; Written by CHANIL PARK (chanil0602@sun.ac.kr) & translate to NCL by EL NOH (noel.modny@gmail.com)
;
; Use      : Solve the quasi-geostrophic (QG) omega equation using the successive over-relaxation method
; Argument : [forcing]  : forcing in the R.H.S of the QG omega equation (3D; nlev x nlat x nlon)
;            [sigma]    : static stability parameter (1D or 3D; nlev or nlev x nlat x nlon)
;            [lon],[lat],[lev] : 1D array corresponding to the longitude, latitude, level which are identical to dimension of [forcing] and [sigma]
;            [ref_lat]  : reference latitude for calculation of f_0. Set as negative for the latitude in the Southern Hemisphere
;            [iternum]  : number of iteration for converged solution
;            [alpha]    : over-relaxation parameter
; Note     : 1) Second-order centered difference is used for spatial differentiation
;            2) Dirichlet boundary condition (omega=0) is applied to the vertical and lateral boundaries
;            3) [lev] should be in the unit of 'Pa'
;            4) Use 1D sigma if you consider area-averaged static stability parameter (only function of [lev]). Otherwise, use 3D sigma varying in the 3D direction (functino of [lon], [lat] and [lev])
;    
; Reference: Park, C., S.-W. Son, and J.-H. Kim, 2021: Role of baroclinic trough in triggering vertical motion during summertime heavy rainfall events in Korea. J. Atmos. Sci.,  78, 1687–1702, https://doi.org/10.1175/JAS-D-20-0216.1.
;            Park, C., S.-W. Son, J. Kim, E.-C. Chang, J.-H. Kim, E. Jo, D.-H. Cha, and S. Jeong, 2021: Diverse synoptic weather patterns of warm-season heavy rainfall events in South Korea. Mon. Wea. Rev., revised.;       
;
;------------------------------------------------------------------------------------------

undef("OmegaEqSolver")
function OmegaEqSolver (forcing:float, sigma:float, lon:float, lat:float, lev:float, ref_lat:float, iternum:integer, alpha:float)
local forcing, sigma, lon, lat, lev, ref_lat, iternum, alpha
begin

; Check of arguments
print("")
if (dimsizes(forcing(:,0,0)).ne.dimsizes(lev) .or. \
    dimsizes(forcing(0,:,0)).ne.dimsizes(lat) .or. \
    dimsizes(forcing(0,0,:)).ne.dimsizes(lon)) then
    print("fatal : Dimension of [forcing] is not identical to [lon], [lat] or [lev]")
end if

if (num(dimsizes(sigma)).eq.1) then
    print("1D sigma is considered...")
if (dimsizes(sigma).ne.dimsizes(lev)) then
    print("fatal : Dimension of [sigma] is not identical to [lev]")
end if
end if

if (num(dimsizes(sigma)).eq.3) then
    print("3D sigma is considered...")
if (dimsizes(sigma(:,0,0)).ne.dimsizes(lev) .or. \
    dimsizes(sigma(0,:,0)).ne.dimsizes(lat) .or. \
    dimsizes(sigma(0,0,:)).ne.dimsizes(lon)) then
    print("fatal : Dimension of [sigma] is not identical to [lev]")
end if
end if

if (ref_lat.lt.-90 .or. ref_lat.gt.90) then
    print("fatal : reference latitude should be set in the range of -90~90")
end if

print("number of iteration: "+iternum)
print("over-relaxation parameter: "+alpha)
print("")

; calculate Coriolis parameter at the reference latitude (f0)
rad   = 4.0*atan(1.0)/180.0 ; radian factor
omega = 7.2921e-5           ; angular velocity of earth
f0    = 2*omega*sin(ref_lat*rad)

; Generalized space differential
nlev = dimsizes(lev)
nlat = dimsizes(lat)
nlon = dimsizes(lon)
dx2  = new((/nlev,nlat,nlon/),float)
dy2  = new((/nlev,nlat,nlon/),float)
dp2  = new((/nlev,nlat,nlon/),float)
dx2  = 0
dy2  = 0
dp2  = 0

dlev = new(nlev,float)
dlat = new(nlat,float)
dlon = new(nlon,float)
dlev = 0
dlat = 0
dlon = 0

dlev(0)        = lev(1)-lev(0)
dlev(1:nlev-2) = (lev(2:)+lev(1:nlev-2))/2 - (lev(0:nlev-3)+lev(1:nlev-2))/2
dlev(nlev-1)   = lev(nlev-1)-lev(nlev-2)

dlat(0)        = lat(1)-lat(0)
dlat(1:nlat-2) = (lat(2:)+lat(1:nlat-2))/2 - (lat(0:nlat-3)+lat(1:nlat-2))/2
dlat(nlat-1)   = lat(nlat-1)-lat(nlat-2)

dlon(0)        = lon(1)-lon(0)
dlon(1:nlon-2) = (lon(2:)+lon(1:nlon-2))/2 - (lon(0:nlon-3)+lon(1:nlon-2))/2
dlon(nlon-1)   = lon(nlon-1)-lon(nlon-2)

dtom = 111000 ; distance on the earth for 1 degree along the great circle (unit:m)
do k = 0,nlev-1
 do j = 0,nlat-1
  dx2(k,j,:) = dlon(:)*cos(lat(j)*rad)*dtom
  dy2(k,j,:) = dlat(j)*dtom
  dp2(k,j,:) = dlev(k)
  ; do i = 0,nlon-1
  ;   dx2(k,j,i) = dlon(i)*cos(lat(j)*rad)*dtom
  ;   dy2(k,j,i) = dlat(j)*dtom
  ;   dp2(k,j,i) = dlev(k)
  ; end do
 end do
end do

; Adjusted space differential
dx = dx2
dy = dy2
dp = dp2
delete(dx2)
delete(dy2)
delete(dp2)
dt2 = new((/nlev,nlat,nlon/),float)
dt2 = 0

;do k = 0,nlev-1
; do j = 0,nlat-1
;  do i = 0,nlon-1
;    if (num(dimsizes(sigma)).eq.1) then
;     dt2(k,j,i) = 2*sigma(k)/(dx(k,j,i)^2) + 2*sigma(k)/(dy(k,j,i)^2) + 2*(f0^2)/(dp(k,j,i)^2)
;    else if (num(dimsizes(sigma)).eq.3) then
;     dt2(k,j,i) = 2*sigma(k,j,i)/(dx(k,j,i)^2) + 2*sigma(k,j,i)/(dy(k,j,i)^2) + 2*(f0^2)/(dp(k,j,i)^2)
;    end if
;    end if
;  end do
; end do
;end do
if (num(dimsizes(sigma)).eq.1) then
 do k = 0,nlev-1
  dt2(k,:,:) = 2*sigma(k)/(dx(k,:,:)^2) + 2*sigma(k)/(dy(k,:,:)^2) + 2*(f0^2)/(dp(k,:,:)^2)
 end do
else if (num(dimsizes(sigma)).eq.3) then
 do k = 0,nlev-1
  dt2(k,:,:) = 2*sigma(k,:,:)/(dx(k,:,:)^2) + 2*sigma(k,:,:)/(dy(k,:,:)^2) + 2*(f0^2)/(dp(k,:,:)^2)
 end do
end if
end if

dt = 1/dt2

; Solve the PDE using the SOR method
W_nume = new(dimsizes(forcing),float) ; array for numerical solution of omega
        ; Dirichlet boundary condition is used (omega=0 for lateral and vertical boundaries)
W_nume = 0

do t = 0,iternum-1
 do k = 1,nlev-2
  do j = 1,nlat-2
   do i = 1,nlon-2
    if (num(dimsizes(sigma)).eq.1) then
     W_nume(k,j,i) = alpha*dt(k,j,i) * ( sigma(k) * (W_nume(k,j,i+1)+(-2)*W_nume(k,j,i)+W_nume(k,j,i-1))/(dx(k,j,i)^2) \
                                       + sigma(k) * (W_nume(k,j+1,i)+(-2)*W_nume(k,j,i)+W_nume(k,j-1,i))/(dy(k,j,i)^2) \
                                       + f0^2     * (W_nume(k+1,j,i)+(-2)*W_nume(k,j,i)+W_nume(k-1,j,i))/(dp(k,j,i)^2) \
                                       - forcing(k,j,i) ) + W_nume(k,j,i)
;     W_nume(k,1:nlat-2,1:nlon-2) = alpha*dt(k,1:nlat-2,1:nlon-2) * ( sigma(k) * (W_nume(k,1:nlat-2,2:nlon-1)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k,1:nlat-2,0:nlon-3))/(dx(k,1:nlat-2,1:nlon-2)^2) \
;                                                                   + sigma(k) * (W_nume(k,2:nlat-1,1:nlon-2)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k,0:nlat-3,1:nlon-2))/(dy(k,1:nlat-2,1:nlon-2)^2) \
;                                                                   + f0^2     * (W_nume(k+1,1:nlat-2,1:nlon-2)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k-1,1:nlat-2,1:nlon-2))/(dp(k,1:nlat-2,1:nlon-2)^2) \
;                                       - forcing(k,1:nlat-2,1:nlon-2) )   + W_nume(k,1:nlat-2,1:nlon-2)
    else if (num(dimsizes(sigma)).eq.3) then
     W_nume(k,j,i) = alpha*dt(k,j,i) * ( sigma(k,j,i) * (W_nume(k,j,i+1)+(-2)*W_nume(k,j,i)+W_nume(k,j,i-1))/(dx(k,j,i)^2) \
                                       + sigma(k,j,i) * (W_nume(k,j+1,i)+(-2)*W_nume(k,j,i)+W_nume(k,j-1,i))/(dy(k,j,i)^2) \
                                       + f0^2         * (W_nume(k+1,j,i)+(-2)*W_nume(k,j,i)+W_nume(k-1,j,i))/(dp(k,j,i)^2) \
                                       - forcing(k,j,i) ) + W_nume(k,j,i)
;     W_nume(k,1:nlat-2,1:nlon-2) = alpha*dt(k,1:nlat-2,1:nlon-2) * ( sigma(k,1:nlat-2,1:nlon-2) * (W_nume(k,1:nlat-2,2:nlon-1)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k,1:nlat-2,0:nlon-3))/(dx(k,1:nlat-2,1:nlon-2)^2) \
;                                                                                 + sigma(k,1:nlat-2,1:nlon-2) * (W_nume(k,2:nlat-1,1:nlon-2)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k,0:nlat-3,1:nlon-2))/(dy(k,1:nlat-2,1:nlon-2)^2) \
;                                                                                 + f0^2            * (W_nume(k+1,1:nlat-2,1:nlon-2)+(-2)*W_nume(k,1:nlat-2,1:nlon-2)+W_nume(k-1,1:nlat-2,1:nlon-2))/(dp(k,1:nlat-2,1:nlon-2)^2) \
;                                       - forcing(k,1:nlat-2,1:nlon-2) )   + W_nume(k,1:nlat-2,1:nlon-2)
    end if
    end if
   end do
  end do
 end do
end do

return (W_nume)

end;/*}}}*/
