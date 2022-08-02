;-- Settings
lev  = 300  ; Pa
type = "x11"
fin1 = "/data3/era5/1.5x1.5/Zpl/Zpl.202112.6hr.grib"
fin2 = "./output/era5_tnflux.nc"

; Load geopotentail height /*{{{*/
undef("load_hgt")
function load_hgt(fnam)
begin
    fin  = addfile(fnam, "r")
    hgt  = fin->Z_GDS0_ISBL
    hgt  = hgt / 9.8    ; geopotentail [m2/s2] -> geopotential height [gpm]
    
    return hgt
end ; /*}}}*/

; Load WAF output (Fx, Fy) /*{{{*/
undef("load_tnflux")
function load_tnflux(fnam)
begin
    f = addfile(fnam,"r")
    Fx = f->Fx
    Fy = f->Fy
    
    return [/ Fx, Fy /]
end ; /*}}}*/

; Plot utility /*{{{*/
undef("plt_waf")
procedure plt_waf(hgt, Fx, Fy, time)
begin
    modulus = sqrt(Fx^2+Fy^2)
    ; lat = hgt&g0_lat_2
    ; lon = hgt&g0_lon_3
    ; modulus = uv2dv_cfd(Fx, Fy, lat, lon, 1)
    copy_VarCoords(Fx, modulus)
    delete_VarAtts(hgt, "long_name")
    delete_VarAtts(hgt, "units")
    delete_VarAtts(Fx, "units")

    wks = gsn_open_wks(type, "waf")
    res = True
    res@gsnFrame = False
    res@gsnDraw = False

    res@tmXBMajorOutwardLengthF = 0.0
    res@tmXBMinorOutwardLengthF = 0.0
    res@tmYLMajorOutwardLengthF = 0.0
    res@tmYLMinorOutwardLengthF = 0.0

    cres = res
    cres@gsnLeftString = "Wave Activity Flux"
    cres@gsnRightString = time
    cres@pmTickMarkDisplayMode = "Always"
    cres@mpGridAndLimbOn       = True  ; Grid on
    cres@mpGridLineDashPattern = 2     ; Grid style
    cres@mpGridLineThicknessF  = 0.5   ; Grid thick
    cres@mpGridLatSpacingF     = 15
    cres@mpGridLonSpacingF     = 15

    cres@cnFillOn = True
    cres@cnLinesOn = False
    cres@cnLevelSelectionMode = "ExplicitLevels"
    cres@cnLevels = (/50, 100, 150, 200, 300, 400/)
    cmap = read_colormap_file("cmocean_matter")
    cmap(0,:) = 0   ; white
    cres@cnFillPalette  = cmap
    ; cres@cnLevels = ispan(-80, 80, 10) * 1.e-5
    ; cres@cnFillPalette = "NCV_blu_red"

    cres@mpLimitMode   = "LatLon"
    cres@mpMinLatF     =  15
    cres@mpMaxLatF     =  80
    cres@mpMinLonF     =  60
    cres@mpMaxLonF     = 160
    cres@pmLabelBarHeightF = 0.06
    cres@pmLabelBarWidthF  = 0.6
    plt = gsn_csm_contour_map(wks, modulus, cres)

    lres = res
    lres@cnLevelSelectionMode = "ExplicitLevels"
    lres@cnLevels = ispan(0, 12000, 60)
    ; lres@cnLevels = ispan(0, 3000, 30)
    lres@cnInfoLabelOn = False
    lres@cnLineLabelPlacementMode = "constant"
    lres@cnLineThicknessF = 2.
    plt1 = gsn_csm_contour(wks, hgt, lres)
    overlay(plt, plt1)

    vres = res
    vres@vcGlyphStyle = "CurlyVector"
    vres@vcMinDistanceF = 0.02
    vres@vcMinMagnitudeF= 10
    vres@vcLineArrowThicknessF = 2.0

    vres@vcRefAnnoPerimOn   = False     ; ref perim off
    vres@vcRefAnnoString2On = False     ; remove the word "reference vector"
    vres@vcRefAnnoParallelPosF   = 1.0  ; NDC coords
    vres@vcRefAnnoOrthogonalPosF = 0.07
    vres@vcRefAnnoFontHeightF    = 0.015
    vres@vcRefLengthF = 0.030
    vres@vcRefMagnitudeF = 100
    plt2 = gsn_csm_vector(wks, Fx, Fy, vres)
    overlay(plt, plt2)
    draw(plt)
    frame(wks)
end ;/*}}}*/

;-- Load data
hgt  = load_hgt(fin1)
era5 = load_tnflux(fin2)
Fx   = era5[0]
Fy   = era5[1]

time = hgt&initial_time0_hours
time:= cd_calendar(time, -3)


;-- Test plot
plt_waf( hgt(0,{lev},:,:), \
         Fx (0,{lev},:,:), \
         Fy (0,{lev},:,:), \
         time(0)           )


