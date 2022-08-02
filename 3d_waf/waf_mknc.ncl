load "./waf_calc.ncl"

;-- 0. Settings
filtering = False
fname = "era5_tnflux.nc"
yyyy  = 2021
mm    = 12
opt = True
opt@plevunit = "hPa"
; opt@masking  = True


; load 6hourly climatology data 1982-2021 /*{{{*/
undef("load_clim")
function load_clim(path, var)
begin
    mon = sprinti("%02d", mm)
    if (var .eq. "uclm") then fin = path+"Upl/"+var+"_1982-2021_"+mon+".nc" end if    
    if (var .eq. "vclm") then fin = path+"Vpl/"+var+"_1982-2021_"+mon+".nc" end if    
    if (var .eq. "tclm") then fin = path+"Tpl/"+var+"_1982-2021_"+mon+".nc" end if    
    if (var .eq. "zclm") then fin = path+"Zpl/"+var+"_1982-2021_"+mon+".nc" end if    
    
    f = addfile(fin,"r")
    data = f->$var$

    if (var .eq. "zclm") then
        data = data / 9.8   ; make geopotential height (gpm)
    end if    
    
    return data
end ; /*}}}*/


; ==================================== Main ==================================== 
;-- 1. Read climatology data
cpath = "/work2/hgkim/2022_gdaps/era5/clim/"
uclm = load_clim(cpath, "uclm") ; [m/s]
vclm = load_clim(cpath, "vclm") ; [m/s]
tclm = load_clim(cpath, "tclm") ; [K]
zclm = load_clim(cpath, "zclm") ; [gpm]

dir  = "/data3/era5/1.5x1.5/Zpl/"
fin  = addfile(dir+"Zpl."+yyyy+""+sprinti("%02d", mm)+".6hr.grib", "r")
hgt  = fin->Z_GDS0_ISBL
hgt  = hgt / 9.8    ; geopotentail [m2/s2] -> geopotential height [gpm]


;-- 2. filtering
; 10-day low-pass fileter
if (filtering) then
    wgt = filwgts_lanczos(20*4+1, 0, 0.1, -999., 1)
    hgt := wgt_runave_n_Wrap(hgt, wgt, 0, 0)
end if    


;-- 3. matching day time 
dims   = dimsizes(hgt)
ntim   = dims(0)
Fx     = new(dims, float)
Fy     = new(dims, float)
Fz     = new(dims, float)
psidev = new(dims, float)
NN     = new(dims, float)

do i = 0, ntim-1
    print(i+" / "+(ntim-1))
    hgt1  =  hgt(i:i,:,:,:)
    uclm1 = uclm(i:i,:,:,:)
    vclm1 = vclm(i:i,:,:,:)
    tclm1 = tclm(i:i,:,:,:)
    zclm1 = zclm(i:i,:,:,:)

    ;-- 4. Calculate TNflux
    tnflux = TNflux(hgt1, uclm1, vclm1, tclm1, zclm1, opt)
    Fx(i,:,:,:)     = tnflux[0]
    Fy(i,:,:,:)     = tnflux[1]
    Fz(i,:,:,:)     = tnflux[2]
    psidev(i,:,:,:) = tnflux[3]
    NN(i,:,:,:)     = tnflux[4]
end do


;-- 5. Save data file
system("rm -f ./output/"+fname)
setfileoption("nc","Format","LargeFile")
fo = addfile("./output/"+fname, "c")
fo->Fx     = Fx    
fo->Fy     = Fy    
fo->Fz     = Fz    
fo->psidev = psidev
fo->NN     = NN    
