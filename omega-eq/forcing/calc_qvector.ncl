load "ncl_myfunction.ncl"

reflat = 38
trunc = True

path = "/work2/hgkim/2022_tc/era5/"
flis1 = systemfunc("ls -1 "+path+"Zpl/Zpl.*.grib")
flis2 = systemfunc("ls -1 "+path+"Tpl/Tpl.*.grib")
g = 9.8 ; m/s2

do i = 0, dimsizes(flis1)-1
    print(i+" / "+(dimsizes(flis1)-1))
    f1   := addfile(flis1(i), "r")
    f2   := addfile(flis2(i), "r")
    Zpl  := f1->Z_GDS0_ISBL;(:,:,{-10:60},{90:190})
    Tpl  := f2->T_GDS0_ISBL;(0,:,{-10:60},{90:190})
    plev := Zpl&lv_ISBL1 * 100    ; Pa
    lat  := Zpl&g0_lat_2
    lon  := Zpl&g0_lon_3
    nlon := dimsizes(lon)

    Tpl  := pot_temp(plev, Tpl, 1, False)

    ;-- 1. RHS: Q-vector
    if (trunc) then
        Tr = 42
        spec_harmonics := shaeC(Zpl)
        spec_harmonics := tri_trunC(spec_harmonics, Tr)
        Zpl = shseC(spec_harmonics, nlon)
        spec_harmonics := shaeC(Tpl)
        spec_harmonics := tri_trunC(spec_harmonics, Tr)
        Tpl := shseC(spec_harmonics, nlon)
    end if    
    
    f  := coriolis_param(lat)
    f0 := coriolis_param(reflat)
    result := qvector1(Zpl, Tpl, f0, plev, 1, 0)
    Qx     := result[0]
    Qy     := result[1]
    Qdiv   := result[2]
    Qdiv_1 := result[5]
    Qdiv_2 := result[8]

    Qdiv   = where(ismissing(Qdiv  ), 0, Qdiv  )
    Qdiv_1 = where(ismissing(Qdiv_1), 0, Qdiv_1)
    Qdiv_2 = where(ismissing(Qdiv_2), 0, Qdiv_2)
    Qdiv@_FillValue   = default_fillvalue(typeof(Qdiv  ))
    Qdiv_1@_FillValue = default_fillvalue(typeof(Qdiv_1))
    Qdiv_2@_FillValue = default_fillvalue(typeof(Qdiv_2))

    Qdiv  (:,0,:,:) = 0
    Qdiv_1(:,0,:,:) = 0
    Qdiv_2(:,0,:,:) = 0

    Qdiv   = where(abs(Qdiv  ).gt.1e-5, 0, Qdiv  )
    Qdiv_1 = where(abs(Qdiv_1).gt.1e-5, 0, Qdiv_1)
    Qdiv_2 = where(abs(Qdiv_2).gt.1e-5, 0, Qdiv_2)
    copy_VarCoords(Zpl, Qdiv)
    copy_VarCoords(Zpl, Qdiv_1)
    copy_VarCoords(Zpl, Qdiv_2)

    opath = path+"omega_solver/qvector/"
    oname = str_split(flis1(i), ".")
    ofile = opath+"Qvector."+oname(1)+"."+oname(2)+".nc"
    system("rm -f "+ofile)
    setfileoption("nc","Format","LargeFile")
    fc = addfile(ofile, "c")
    fc->Qx = Qx
    fc->Qy = Qy
    fc->Qfor    = Qdiv
    fc->Qparfor = Qdiv_1
    fc->Qperfor = Qdiv_2
end do    

