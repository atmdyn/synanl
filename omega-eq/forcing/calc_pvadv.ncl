load "ncl_myfunction.ncl"

reflat = 38
trunc = True

path = "/work2/hgkim/2022_tc/era5/"
flis1 = systemfunc("ls -1 "+path+"Zpl/Zpl.*.grib")
g = 9.8 ; m/s2
rad = 4.0*atan(1.0)/180.0 ; radian factor
omega=7.2921e-5 ; angular veloity of earth
a=6.37122e6     ; raidus of earth

do i = 0, dimsizes(flis1)-1
    print(i+" / "+(dimsizes(flis1)-1))
    f1   := addfile(flis1(i), "r")
    Zpl  := f1->Z_GDS0_ISBL;(:,:,{-10:60},{90:190})
    plev := Zpl&lv_ISBL1 * 100    ; Pa
    lat  := Zpl&g0_lat_2
    lon  := Zpl&g0_lon_3
    nlon := dimsizes(lon)

    ;-- 2. RHS: planetary vorticity advection by thermal wind
    if (trunc) then
        Tr = 42
        spec_harmonics := shaeC(Zpl)
        spec_harmonics := tri_trunC(spec_harmonics, Tr)
        Zpl = shseC(spec_harmonics, nlon)
    end if    
    

    ; beta (df/dy)
    f  = 2*omega*sin(lat*rad) ; Coriolis parameter
    f0 = 2*omega*sin(reflat*rad) ; Coriolis parameter for reference latitude
    dy = a*abs(lat(1)-lat(0))*rad
    beta = center_finite_diff(f,dy,False,0)
    beta0 = beta(minind(lat-reflat))

    ; meridional geostrophic wind
    vg := new(dimsizes(Zpl),"float")
    dx := new((/dimsizes(lat)/),"float")
    do j=0,dimsizes(lat)-1
     dx(j) = a*abs(cos(lat(j)*rad)) * abs(lon(1)-lon(0))*rad
    end do
    do j=0,dimsizes(lat)-1
      vg(:,:,j:j,:) = 1/f0 * center_finite_diff_n(Zpl(:,:,j:j,:),dx(j),True,0,3)
    end do

    ; vertical derivative of the meridional geostrophic wind
    dvgdp := center_finite_diff_n(vg,plev,False,2,1)
    BETAeffect := f0 * beta0 * dvgdp
    copy_VarCoords(Zpl,BETAeffect)
    BETAeffect = where(ismissing(BETAeffect),0,BETAeffect)
    BETAeffect@_FillValue = -999.

    BETAeffect = where(abs(BETAeffect).gt.1e-5,0,BETAeffect)

    opath = path+"omega_solver/pvadv/"
    oname = str_split(flis1(i), ".")
    ofile = opath+"pvadv."+oname(1)+"."+oname(2)+".nc"
    system("rm -f "+ofile)
    setfileoption("nc","Format","LargeFile")
    fc = addfile(ofile, "c")
    fc->betafor = BETAeffect
end do    

