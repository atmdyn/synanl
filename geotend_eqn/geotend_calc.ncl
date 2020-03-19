load "./src_test/poisson_solver.ncl"
load "./geotend_func.ncl"

dtype = "fms"   ; era(ERA-Interim), fms(FMS)

; Read file ==================================================/*{{{*/
;-- Era-Interim
if (dtype .eq. "era") then
    y4m2 = "201602"
    path = "/data2/eraint/1.5x1.5/"
    fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
    U  = fU->U_GDS0_ISBL 
    fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
    V  = fV->V_GDS0_ISBL 
    fin  = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
    T    = fin->T_GDS0_ISBL(0:11,{200:1000},:,:)    ; Temperature [K]
    fin  = addfile(path+"Zpl/"+"Zpl."+y4m2+".6hr.grib", "r")
    PHI  = fin->Z_GDS0_ISBL(0:11,{200:1000},:,:)    ; Geopotential [m**2 s**-2]
    fin  = addfile(path+"PHYpl/"+"Tphy.pl."+y4m2+".nc", "r")
    J_cp = fin->Tphy(0:11,{200:1000},:,:)/(3600.*6) ; =J/Cp [K/6hour] -> [K/s]
else
    ;-- FMS output
    path = "../data/"
    ; fin  = addfile(path+"fms_HS_U.nc", "r")
    ; U    = fin->U(0:4,{200:1000},:,:)
    ; fin  = addfile(path+"fms_HS_V.nc", "r")
    ; V    = fin->V(0:4,{200:1000},:,:)
    fin  = addfile(path+"fms_HS_T.nc", "r")
    T    = fin->T(0:9,{200:1000},:,:)       ; Temperature [K]
    fin  = addfile(path+"fms_HS_Z.nc", "r")
    PHI  = fin->Z(0:9,{200:1000},:,:)*9.8   ; Geopotential [m**2 s**-2]
    copy_VarCoords(T, PHI)
    fin  = addfile(path+"fms_HS_Tdt.nc", "r")
    J_cp = fin->Tdt(0:9,{200:1000},:,:)     ; [K/s]
end if    

; J_cp = J_cp*86400.
; wks = gsn_open_wks("x11", "nonamed")
; res = True
; res@cnFillOn = True
; plt = gsn_csm_contour_map(wks, J_cp(1,{500},:,:), res)
; exit
;============================================================/*}}}*/


dim  = dimsizes(PHI)
; PHId = dim_avg_n(reshape(PHI, (/3, 4, dim(1),dim(2),dim(3)/)),1)
; X = center_finite_diff_n(PHId, 86400., False, 1, 0)
; X!0 = "days"
; copy_VarCoords(T(0,:,:,:), X(0,:,:,:))
X = center_finite_diff_n(PHI, 86400., False, 1, 0)
; X = center_finite_diff_n(PHI, 21600., False, 1, 0)
copy_VarCoords(T, X)

F  = new(dim, "float")
F1 = new(dim, "float")
F2 = new(dim, "float")
F3 = new(dim, "float")
copy_VarCoords(T, F)
copy_VarCoords(T, F1)
copy_VarCoords(T, F2)
copy_VarCoords(T, F3)
solution = new(dim, "double")
copy_VarCoords(T, solution)

do tidx = 0, 9
    print("tidx : "+tidx)
    opt1 = True
    opt1@T = T(tidx,:,:,:)
    dummy = geotend_forcing(PHI(tidx,:,:,:), J_cp(tidx,:,:,:), opt1)
    F1(tidx,:,:,:) = dummy@F1
    F2(tidx,:,:,:) = dummy@F2
    F3(tidx,:,:,:) = dummy@F3
    delete_VarAtts(dummy, "F1")
    delete_VarAtts(dummy, "F2")
    delete_VarAtts(dummy, "F3")
    F(tidx,:,:,:) = dummy

    ; Solve equation using SOR(Successive Over Relexation) method
    opt = True
    opt@T0 = dim_avg_Wrap(dim_avg_Wrap(T(tidx,:,:,:)))
    opt@lev_units = "hPa"
    opt@dt = 1.e8
    opt@f0 = 1.e-4
    opt@thres = 1.e-10
    opt@upper_bound = 0.
    opt@lower_bound = 0.
    solution(tidx,:,:,:) = poisson_solver(F(tidx,:,:,:), opt)
end do


ofile = "./geotend.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->F1  = F1
fc->F2  = F2
fc->F3  = F3
fc->X   = X
fc->solution = solution
fc->forcing = F
















