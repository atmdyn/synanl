load "ncl_myfunction.ncl"

reflat = 38
trunc  = True

path  = "/work2/hgkim/2022_tc/era5/"
flis1 = systemfunc("ls -1 "+path+"Tphy/Tphy.*.nc")
g     = 9.8 ; m/s2
rad   = 4.0*atan(1.0)/180.0 ; radian factor
omega = 7.2921e-5           ; angular veloity of earth
a     = 6.37122e6           ; raidus of earth
R     = 287.                ; Gas constant
Cp    = 1005.               ; Specific heat at constant pressure [J/kg/K]

do i = 0, dimsizes(flis1)-1
    print(i+" / "+(dimsizes(flis1)-1))
    f1   := addfile(flis1(i), "r")
    Tphy := f1->Q1;(:,:,{-10:60},{90:190})
    plev := Tphy&lv_ISBL1 * 100    ; Pa
    lat  := Tphy&g0_lat_2
    lon  := Tphy&g0_lon_3
    nlon := dimsizes(lon)

    ;-- 2. RHS: planetary vorticity advection by thermal wind
    if (trunc) then
        Tr = 128
        Tphy = where(ismissing(Tphy),0,Tphy)
        spec_harmonics := shaeC(Tphy)
        spec_harmonics := tri_trunC(spec_harmonics, Tr)
        Tphy = shseC(spec_harmonics, nlon)
    end if    
    
    
    dum := tofloat(Tphy/Cp) ; J/kg/s -> K/s
    copy_VarCoords(Tphy, dum)
    Tphy := dum
    copy_VarCoords(dum, Tphy)

    ; Laplacian of total physics tendency (K/s)
    dx = new((/dimsizes(lat)/),"float")
    do j=0,dimsizes(lat)-1
     dx(j) = a*abs(cos(lat(j)*rad)) * abs(lon(1)-lon(0))*rad
    end do
    dy = a*abs(lat(1)-lat(0))*rad

    d2Tphydx2 := Tphy*0
    do j=0,dimsizes(lat)-1
     d2Tphydx2(:,:,j:j,:) = center_finite_diff_n(Tphy(:,:,j:j,:),dx(j),False,2,3)
    end do
    do j=0,dimsizes(lat)-1
     d2Tphydx2(:,:,j:j,:) = center_finite_diff_n(d2Tphydx2(:,:,j:j,:),dx(j),False,2,3)
    end do

    d2Tphydy2 := Tphy*0
    d2Tphydy2 = center_finite_diff_n(Tphy,dy,False,2,2)
    d2Tphydy2 = center_finite_diff_n(d2Tphydy2,dy,False,2,2)

    LapTphy := Tphy*0
    LapTphy = (d2Tphydx2 + d2Tphydy2)

    ; final form of diabatic forcing term
    do k=0,dimsizes(plev)-1
     LapTphy(:,k,:,:) = (-1./plev(k)) * (R/Cp) * LapTphy(:,k,:,:)*Cp
    end do
    copy_VarCoords(Tphy,LapTphy)
    LapTphy = where(ismissing(LapTphy),0,LapTphy)
    LapTphy@_FillValue = -999.

    LapTphy = where(abs(LapTphy).gt.1e-5,0,LapTphy)

    opath = path+"omega_solver/diab/"
    oname = str_split(flis1(i), ".")
    ofile = opath+"diab."+oname(1)+"."+oname(2)+".nc"
    system("rm -f "+ofile)
    setfileoption("nc","Format","LargeFile")
    fc = addfile(ofile, "c")
    fc->diafor = LapTphy

end do    

