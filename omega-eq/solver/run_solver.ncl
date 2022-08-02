load "omega_solver.ncl"
load "../../2022_tc.ncl"

; calc_refatm_sigma /*{{{*/
undef("calc_refatm_sigma")
function calc_refatm_sigma(p, t0, pdim)
local R, ss0, sigma
begin
    R = 287     ; dry gas constant

    ; Compute sigma for reference atmosphere (See Holton 5th 209p)
    sigma = conform(t0, R/p, pdim)*static_stability(p, t0, pdim, 0)
    copy_VarCoords(t0, sigma)
    sigma = where(ismissing(sigma), 0, sigma)
    sigma@_FillValue = default_fillvalue(typeof(sigma))

    return sigma
end ; /*}}}*/


som_y = 1
som_x = 3
baselat = 25
reflat = 38
alpha   = 1.5    ; over-relaxation parameter
iternum = 500    ; specified number of iteration
cluster_case = 3 ; 1, 2, 3


;-- data
name = call_tc_list(1985, 2018)
name := get_JJASO_tc(name)
name := get_som_tc(name, som_y, som_x, cluster_case)
t0   = get_t0_tc(name, baselat)
qsort(t0)

;-- 500 hPa synop data
yr = cd_calendar(t0, -1)
yr:= tostring( get_unique_values(yr) )
path = "/work2/hgkim/2022_tc/era5/"
flis1 = path+"omega_solver/qvector/Qvector."+yr+".6hr.nc"
flis2 = path+"omega_solver/pvadv/pvadv."+yr+".6hr.nc"
flis3 = path+"omega_solver/diab/diab."+yr+".6hr.nc"
flis4 = path+"Tpl/Tpl."+yr+".6hr.grib"
flis5 = path+"Wpl/Wpl."+yr+".6hr.grib"


f    := addfiles(flis4, "r")
temp := f[:]->T_GDS0_ISBL({t0},:,{-10:60},{90:190})  ; K
temp := dim_avg_n_Wrap(temp, 0)
lon  := temp&g0_lon_3
lat  := temp&g0_lat_2
prs  := temp&lv_ISBL1*100.     ; hPa -> Pa

; Compute sigma for reference atmosphere (See Holton 5th 209p)
temp0 := area_average(temp)
sigma := calc_refatm_sigma(prs, temp0, 0)

; Q-vector forcing
f1   := addfiles(flis1, "r")
Qfor := f1[:]->Qfor({t0},:,{-10:60},{90:190})
Qfor := dim_avg_n_Wrap(Qfor, 0)
Qvec1:= dim_avg_n_Wrap(f1[:]->Qparfor({t0},:,{-10:60},{90:190}), 0)
Qvec2:= dim_avg_n_Wrap(f1[:]->Qperfor({t0},:,{-10:60},{90:190}), 0)
Qx   := dim_avg_n_Wrap(f1[:]->Qx({t0},:,{-10:60},{90:190}), 0)
Qy   := dim_avg_n_Wrap(f1[:]->Qy({t0},:,{-10:60},{90:190}), 0)

; beta forcing
f2   := addfiles(flis2, "r")
bfor := f2[:]->betafor({t0},:,{-10:60},{90:190})
bfor := dim_avg_n_Wrap(bfor, 0)

; diab forcing
f3   := addfiles(flis3, "r")
dfor := f3[:]->diafor({t0},:,{-10:60},{90:190})
dfor := dim_avg_n_Wrap(dfor, 0)
dfor := tofloat(dfor)
copy_VarCoords(bfor, dfor)

; solve omega
fors   = Qfor + bfor + dfor
W      = OmegaEqSolver(fors, sigma, lon, lat, prs, reflat, iternum, alpha)
W_Qvec = OmegaEqSolver(Qfor, sigma, lon, lat, prs, reflat, iternum, alpha)
W_beta = OmegaEqSolver(bfor, sigma, lon, lat, prs, reflat, iternum, alpha)
W_diab = OmegaEqSolver(dfor, sigma, lon, lat, prs, reflat, iternum, alpha)
copy_VarCoords(temp, W)
copy_VarCoords(temp, W_Qvec)
copy_VarCoords(temp, W_beta)
copy_VarCoords(temp, W_diab)

f5   := addfiles(flis5, "r")
omg   = f5[:]->W_GDS0_ISBL({t0},:,{-10:60},{90:190})
omg  := dim_avg_n_Wrap(omg, 0)

opath = path+"omega_solver/omega/"
ofile = opath+"omega."+som_y+"x"+som_x+"_"+cluster_case+".nc"
system("rm -f "+ofile)
setfileoption("nc","Format","LargeFile")
fc = addfile(ofile, "c")
fc->W_all  = W
fc->W_Qvec = W_Qvec
fc->W_beta = W_beta
fc->W_diab = W_diab
fc->sigma = sigma
fc->Qx   = Qx
fc->Qy   = Qy
fc->Qfor = Qfor
fc->Qparfor= Qvec1
fc->Qperfor= Qvec2
fc->bfor = bfor
fc->dfor = dfor
fc->W_org = omg
