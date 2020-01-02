; ------------------------------------------------ ;
; Created on  Thu Jan 02  18:37:01 2020
; @author : Hyeong-Gyu Kim.
;
;-- Description
;  - Full version of apparent_calc.ncl and apparent_mknc.ncl
;  - For refering.
; 
;-- Update
;  - 2020-01-02 : Maked by ADLab
;
;-- Reference
;  - Yanai et al., (1973)
;  - Lee et al., (2008)
;  - Fueglistaler et al., (2009)
; ------------------------------------------------ ;

y4m2 = "201607"
path = "/data2/eraint/1.5x1.5/"

; Read file
fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
U  = fU->U_GDS0_ISBL 
fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
V  = fV->V_GDS0_ISBL 
fW = addfile(path+"Wpl/"+"Wpl."+y4m2+".6hr.grib", "r")
W  = fW->W_GDS0_ISBL 
fT = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T  = fT->T_GDS0_ISBL 
prs = T&lv_ISBL1*100. ; pressure [Pa]
lat = T&g0_lat_2 ; latitude
ny = dimsizes(lat)

fQ = addfile(path+"Qpl/"+"Qpl."+y4m2+".6hr.grib", "r")
Q  = fQ->Q_GDS0_ISBL

; Calc equation
a  = 6.37122e06         ; radius of the earth (NCL function)
R  = 287.               ; [J/kg/K]
cp = 1005.              ; [J/kg/K]
L  = 2.501e6            ; [J/kg]
; k  = 2./7.  ; =R/cp
k  = R/cp
dt = 6*3600.            ; [s]
pi  = 4*atan(1.)        ; = 3.141592
d2r = pi/180.
; dx  = 1.5*d2r*a*cos(lat) ; 1.5 degree
dTdt = center_finite_diff_n(T, dt, False, 1, 0)
dTdx = dTdt
do iy = 0, ny-1
    dx = 1.5*d2r*a*cos(lat(iy)*d2r) ; 1.5 degree
    dTdx(:,:,iy,:) = center_finite_diff_n(T(:,:,iy,:), dx, True, 1, 2)
end do
dy   = lat*d2r*a          ; 1.5 degree
dTdy = center_finite_diff_n(T, dy, False, 1, 2)

; Finite difference in log-P
; d/dp => 1/p*d/dlogp. otherwise just use p for vert. diff.
logp = log(prs)
dTdp = center_finite_diff_n(T, logp, False, 1, 1)
dTdp = dTdp/conform(dTdp, prs, 1)          ; take it back to pressure

; Advection term
Tadv_x = U*dTdx
Tadv_y = V*dTdy
Tadv = Tadv_x + Tadv_y

; Vertical advection (THETA)
Tadv_p = -W*(k*T/conform(T, prs, 1) - dTdp)

Q1 = cp*(dTdt + Tadv + Tadv_p)
copy_VarCoords(T, Q1)
printVarSummary(Q1)
Q1@long_name = "Apparent Heat Source"
Q1@units = "J/kg/s"

; Q2 part
dqdt = center_finite_diff_n(Q, dt, False, 1, 0)
dqdx = dqdt
do iy = 0, ny-1
    dx = 1.5*d2r*a*cos(lat(iy)*d2r) ; 1.5 degree
    dqdx(:,:,iy,:) = center_finite_diff_n(Q(:,:,iy,:), dx, True, 1, 2)
end do
dqdy = center_finite_diff_n(Q, dy, True, 1, 2)
dqdp = center_finite_diff_n(Q, logp, False, 1, 1)
dqdp = dqdp/conform(dqdp, prs, 1)          ; take it back to pressure

qadv_x = U*dqdx
qadv_y = V*dqdy
qadv   = qadv_x + qadv_y
qadv_p = W*dqdp

Q2 = -L*(dqdt + qadv + qadv_p)
copy_VarCoords(T, Q2)
printVarSummary(Q2)
Q2@long_name = "Apparent Moisture Sink"
Q2@units = "J/kg/s"

setfileoption("nc","Format","LargeFile")
ofile = "apparent_heat_moist.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->Q1=Q1
fc->Q2=Q2
fc->Tadv_x = Tadv_x
fc->Tadv_y = Tadv_y
fc->Tadv_p = Tadv_p
fc->dTdt   = dTdt
