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

path = "/work2/hgkim/2022_tc/era5/Tpl/"
flis = systemfunc("ls -1 "+path+"Tpl.*.grib")

R = 287     ; dry gas constant
do i = 0, dimsizes(flis)-1
    f := addfile(flis(i), "r")
    temp := f->T_GDS0_ISBL(:,:,{-10:60},{90:190})  ; K
    prs  := f->lv_ISBL1*100.     ; hPa -> Pa

    ; Compute sigma for reference atmosphere (See Holton 5th 209p)
    t0  = area_average(temp)
    sigma = calc_refatm_sigma(prs, t0, 1)
end do    

