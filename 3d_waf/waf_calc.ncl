; ERA5  /*{{{*/
undef("load_era5")
function load_era5(var, lev, fcst)
begin

    ddir = "/data3/era5/1.5x1.5/"+var+"/"
    flis = var + (/".202007.6hr.grib", ".202008.6hr.grib", ".202009.6hr.grib"/)
    f = addfiles(ddir+flis, "r")
    if (var .eq. "Upl") then fvar = "U_GDS0_ISBL" end if    
    if (var .eq. "Vpl") then fvar = "V_GDS0_ISBL" end if    
    if (var .eq. "Zpl") then fvar = "Z_GDS0_ISBL" end if    
    if (var .eq. "Tpl") then fvar = "T_GDS0_ISBL" end if    
    if (var .eq. "Qpl") then fvar = "Q_GDS0_ISBL" end if    
    
    units = "hours since 1800-01-01 00:00"
    stim  = cd_inv_calendar(2020, 07, 27,  0, 0, 0, units, 0) + fcst
    etim  = cd_inv_calendar(2020, 08, 31, 12, 0, 0, units, 0) + fcst

    if (lev .eq. -1) then
        kimlev = (/ 1000, 975, 950, 925, 900, 875, 850, 800, 750, \
                     700, 650, 600, 550, 500, 450, 400, 350, 300, \
                     250, 200, 150, 100,  70,  50,  30,  20,  10, \
                       7,   5,   3,   1 /)
        data = f[:]->$fvar$({stim:etim:2},{kimlev},:,:)
    else
        data = f[:]->$fvar$({stim:etim:2},{lev},:,:)
    end if    

    if (var .eq. "Zpl") then 
        data = data/9.8 
        data@long_name = "Geopotential height"
        data@units = "m"
    end if    ; m2/s2 -> m
    if (var .eq. "Qpl") then 
        data = data*1e3 
        data@units = "g/kg"
    end if    ; kg/kg -> g/kg

    return data
end ; /*}}}*/

; get_dims /*{{{*/
undef("get_dims")
function get_dims(var:numeric)
; 1. Description
;   - Assume `var` dimension is [...] x [time] x [lev] x [lat] x [lon]
; 
; 2. Input
;   - Variables >= 1 dimension (with coordinate info)
; 
; 3. Output
;   - dimension size of variables
;
; 4. history
;   2021-10-18 : Make & test function (by Hyeong-Gyu Kim)
begin

    ; check dimension
    dims = dimsizes(dimsizes(var))
    ndim = dimsizes(var)

    if (dims .ge. 4) then 
        tim_dim_idx   = dims - 4
        dims@nt       = ndim(tim_dim_idx)
        if ( isdimnamed(var, tim_dim_idx) ) then
            tim_name  = var!tim_dim_idx
            dims@tim_name = tim_name
            if ( iscoord(var, var!tim_dim_idx) ) then
                dims@tim = var&$tim_name$
            end if    
        end if    
    end if

    if (dims .ge. 3) then
        lev_dim_idx   = dims - 3
        dims@np       = ndim(lev_dim_idx)
        if ( isdimnamed(var, lev_dim_idx) ) then
            lev_name  = var!lev_dim_idx
            dims@lev_name = lev_name
            if ( iscoord(var, var!lev_dim_idx) ) then
                dims@lev = var&$lev_name$
            end if    
        end if    
    end if    

    if (dims .ge. 2) then
        lat_dim_idx   = dims - 2 
        dims@ny       = ndim(lat_dim_idx)
        if ( isdimnamed(var, lat_dim_idx) ) then
            lat_name  = var!lat_dim_idx
            dims@lat_name = lat_name
            if ( iscoord(var, var!lat_dim_idx) ) then
                dims@lat = var&$lat_name$
            end if    
        end if    
    end if    

    if (dims .ge. 1) then
        lon_dim_idx   = dims - 1
        dims@nx       = ndim(lon_dim_idx)
        if ( isdimnamed(var, lon_dim_idx) ) then
            lon_name  = var!lon_dim_idx
            dims@lon_name = lon_name
            if ( iscoord(var, var!lon_dim_idx) ) then
                dims@lon = var&$lon_name$
            end if    
        end if    
    end if    
    
    return dims
end ; /*}}}*/

; Calculate Takaya & Nakamura (2001) 3D Wave Activity flux  /*{{{*/
undef("TNflux")
function TNflux(  hgt[*][*][*][*]:numeric, \ ; geopotential height [m]
                 uclm[*][*][*][*]:numeric, \ ; climatology of zonal wind [m/s]
                 vclm[*][*][*][*]:numeric, \ ; climatology of meridional wind [m/s]
                 tclm[*][*][*][*]:numeric, \ ; climatology of air temperature [deg K]
                 zclm[*][*][*][*]:numeric, \ ; climatology of geopotential height [m]
                           opt:logical )
; 1. Description
; 3-D wave-activity flux derived by Takaya and Nakamura (1999, 2001)
; See (38) of Takaya and Nakamura (2001, JAS)
;
; - Basic state  : daily(or hourly) climatology
; - Perturbation : deviation from climatology for each calendar day
;
; Calculating the flux for each time of input data.
;
; 2. Arguments
; - hgt  : geopotential height [m]                
; - uclm : climatology of zonal wind [m/s]
; - vclm : climatology of meridional wind [m/s]
; - tclm : climatology of air temperature [deg K]
; - zclm : climatology of geopotential height [m]
; -  opt : opt@plevunit = "hPa" : (default)
; -  opt : opt@plevunit = "Pa"  : The unit of level change [Pa] -> [hPa]
;          opt@masking  = False : no masking (default)
;          opt@masking  = True  : masking when TNflux < 5 [m/s]
;
; 3. Output
; - TN flux of x, y, z-component, [m2/s2]
; - psidev: daily QG stream function anomaly [m2/s]
; - NN: Brunt Vaisala frequency [1/s2]
;
; Note!  
; (i)  10-day low-pass filter is applied
;       (Assume that phase speed (Cp) is almost zero)
; (ii) "p" (pressure) in Eq. (5.7) is dimensionless (pressure/1000mb)
; (iii) The year 2000 is a leap year (366 days)
; (iv)  Using height data of 1999 and 2001 is for time filtering. 
;       The filter makes the first and last several days unavailable.
begin
    ; Option
    plevunit = "hPa"
    masking  = False
    if (isatt(opt, "plevunit")) then plevunit = opt@plevunit end if
    if (isatt(opt,  "masking")) then  masking = opt@masking  end if

    ; Set parameter
    gc = 290           ; Gas constant
    ga = 9.80665       ; Gravitational acceleration [m/s2]
    re = 6378388       ; Radius of the earth [m]
    sclhgt = 8000.     ; scale height [m]
    pi = atan(1.0)*4.  ; = 3.141592

    ; Get coordinate info
    dims = get_dims(hgt)
    ; time = dims@tim
    lev  = dims@lev
    lat  = dims@lat
    lon  = dims@lon
    dim  = dimsizes(hgt)
    nt   = dim(0)
    ; nt   = dimsizes(time)
    np   = dimsizes(lev)
    ny   = dimsizes(lat)
    nx   = dimsizes(lon)

    if (plevunit .eq. "Pa") then lev := lev / 100. end if    

    ; height anomaly
    zano = hgt - zclm 

    ; Coriolis parameter
    f = 2.*2.*pi/(60.*60.*24.)*sin(pi/180. * lat(:))
    f!0 = "lat"
    f&lat = lat
    f@_FillValue = default_fillvalue(typeof(f))
    f = where( abs(lat) .lt. 10., f@_FillValue, f) ; missing for 10S - 10N

    ; cosine
    clat = cos( lat(:)*pi/180. )

    ; 1-D -> 4-D
    lev4d  = conform_dims( dimsizes(zano),  lev, 1 )
    clat4d = conform_dims( dimsizes(zano), clat, 2 )
    f4d    = conform_dims( dimsizes(zano),    f, 2 )

    ; vertical gradient of potential temperature [K/m]
    dthetadz = center_finite_diff_n( tclm*(1000./lev4d )^0.286, \
                                         -sclhgt*log(lev/1000.), \
                                                    False, 0, 1 )

    ; Brunt Vaisala frequency
    NN = (gc*(lev4d /1000.)^0.286)/sclhgt * dthetadz
    NN@_FillValue = default_fillvalue(typeof(NN))
    NN = where(NN .gt. 0, NN, NN@_FillValue)

    ; magnitude of climatological wind
    cumag = sqrt(uclm^2 + vclm^2)
    cumag@_FillValue = default_fillvalue(typeof(cumag))
    cumag = where(cumag .gt. 0, cumag, cumag@_FillValue)

    ; QG steam function for anomaly
    psidev = zano*ga / f4d 

    ;dpsidev/dlon
    dpsidevdlon = center_finite_diff_n( psidev, \
                                   lon*pi/180., \
                                     True, 0, 3 )

    ;ddpsidev/dlonlon
    ddpsidevdlonlon = center_finite_diff_n( dpsidevdlon, \
                                            lon*pi/180., \
                                              True, 0, 3 )

    ;dpsidev/dlat
    dpsidevdlat = center_finite_diff_n( psidev, \
                                   lat*pi/180., \
                                    False, 0, 2 )

    ;ddpsidev/dlonlat
    ddpsidevdlonlat = center_finite_diff_n( dpsidevdlon, \
                                            lat*pi/180., \
                                             False, 0, 2 )

    ;ddpsidev/dlatdlat
    ddpsidevdlatlat = center_finite_diff_n( dpsidevdlat, \
                                            lat*pi/180., \
                                             False, 0, 2 )

    ;dpsidev/dz
    dpsidevdz = center_finite_diff_n( psidev, \
                       -sclhgt*log(lev/1000.), \
                                  False, 0, 1 )

    ;ddpsidev/dlondz
    ddpsidevdlonz = center_finite_diff_n( dpsidevdlon, \
                                -sclhgt*log(lev/1000.), \
                                           False, 0, 1 )

    ;ddpsidev/dlatdz
    ddpsidevdlatz = center_finite_diff_n( dpsidevdlat, \
                                -sclhgt*log(lev/1000.), \
                                           False, 0, 1 )

    xuterm = (dpsidevdlon*dpsidevdlon - psidev*ddpsidevdlonlon)
    xvterm = (dpsidevdlon*dpsidevdlat - psidev*ddpsidevdlonlat)
    ;yuterm = xvterm
    yvterm = (dpsidevdlat*dpsidevdlat - psidev*ddpsidevdlatlat)

    zuterm = (dpsidevdlon*dpsidevdz - psidev*ddpsidevdlonz)
    zvterm = (dpsidevdlat*dpsidevdz - psidev*ddpsidevdlatz)

    ; x-component of (38)
    Fx = lev4d /1000./(2.*cumag*re*re)*( uclm/clat4d  * xuterm + vclm * xvterm )

    ; y-component 
    Fy = lev4d /1000./(2.*cumag*re*re)*( uclm*xvterm + clat4d *vclm*yvterm )

    ; z-component 
    Fz = lev4d /1000.*f4d *f4d /(2.*cumag*NN*re)*(uclm*zuterm + clat4d *vclm*zvterm)

    ; Mask out where westerlies is small or negative (less than 5 m/s). by using mask
    if (masking) then
        Fx = mask( Fx, uclm.lt.5, False )
        Fy = mask( Fy, uclm.lt.5, False )
        Fz = mask( Fz, uclm.lt.5, False )
    end if

    ; For output
    Fx     := tofloat(Fx    )
    Fy     := tofloat(Fy    )
    Fz     := tofloat(Fz    )
    psidev := tofloat(psidev)
    NN     := tofloat(NN    )
    copy_VarCoords(hgt, Fx)
    copy_VarCoords(hgt, Fy)
    copy_VarCoords(hgt, Fz)
    copy_VarCoords(hgt, psidev)
    copy_VarCoords(hgt, NN)

    Fx@units     = "m^2/s^2"
    Fy@units     = "m^2/s^2"
    Fz@units     = "m^2/s^2"
    psidev@units = "m^2/s"
    NN@units     = "1/s^2"

    return [/     Fx, \ ; Wave activity x-component [m2/s2]
                  Fy, \ ; Wave activity y-component [m2/s2]
                  Fz, \ ; Wave activity z-component [m2/s2]
              psidev, \ ; daily QG stream function anomaly [m2/s]
                  NN /] ; NN: Brunt Vaisala frequency [1/s2]
end ; /*}}}*/
