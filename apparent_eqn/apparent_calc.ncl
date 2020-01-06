; ------------------------------------------------ ;
; Created on  Thu Jan 02  18:37:01 2020
; @author : Joowan Kim, Hyeong-Gyu Kim. Hyeon-Seok Do
;
;-- Description
;  - Function code for computing Q1(apparent heat source) and Q2(apparent moisture sink)
;  - Input variable must be 4 dimension([Time]x[plev]x[lat]x[lon])
;  - In this code, [Pa] of pressure unit is used.
;  - Note! You must set opt@prs_unit
;    ex) opt@prs_unit = "hPa"
;        opt@prs_unit = "Pa"
;
;
;-- Nomenclature
;  - U, V     : zonal, meridional wind components [m/s]
;  - q        : specific humidity                 [g/kg]
;  - T        : temperature                       [K] or [C]
;  - W(omega) : vertical velocity                 [Pa/s]
;  - opt      : set to True if pressure units is "hPa" (opt@prs_unit = "hPa)
; 
;
;-- Output
;  - Q1 : Apparent heat source   [J/kg/s]
;  - Q2 : Apparent moisture sink [J/kg/s]
;   Q1 = Cp*(dTdt + Tadv + Tadv_p)
;   Q2 = -L*(dqdt + qadv + qadv_p)
;
;
;-- Update
;  - 2020-01-02 : Maked by ADLab
;
;
;-- Reference
;  - Yanai et al., 1973:
;   Determination of bulk properties of tropical cloud clusters from 
;   large-scale heat and moisture budgets. Journal of the Atmospheric 
;   Sciences, 30(4), 611-627.
;   https://doi.org/10.1175/1520-0469(1973)030<0611:DOBPOT>2.0.CO;2
;
;  - Lee et al., (2008)
;   Lee et al., 2009:
;   Heavy Rainfall Events Lasting 18 Days from July 31 to August 17, 
;   1998, over Korea. J. Meteorol. Soc. Japan, 86, 313–333.
;   https://doi.org/10.2151/jmsj.86.313
;
;  - Fueglistaler et al., 2009:
;   The diabatic heat budget of the upper troposphere and lower/mid 
;   stratosphere in ECMWF reanalyses. Q. J. R. Meteorol. Soc., 135, 21–37
;   https://doi.org/10.1002/qj.361
;
;
;-- Further Reading:
;  - https://renqlsysu.github.io/2019/02/01/apparent_heat_source/
;
;  - Yanai, M. (1961): 
;   A detailed analysis of typhoon formation.
;   J. Meteor. Soc. Japan , 39 , 187��214
;   https://www.jstage.jst.go.jp/article/jmsj1923/39/4/39_4_187/_pdf
;
;  - Yanai, M and T.Tomita (1998):
;   https://pdfs.semanticscholar.org/fb57/a6a59cc4a684194b5e622ea6f875d0b4439a.pdf
; ------------------------------------------------ ;

; Calculate apparent heat source (Q1) /*{{{*/
undef("apparent_heat_source") 
function apparent_heat_source(T[*][*][*][*]:numeric, \  ; Temperature     [K]
                              U[*][*][*][*]:numeric, \  ; Zonal wind      [m/s]
                              V[*][*][*][*]:numeric, \  ; Meridional wind [m/s]
                              W[*][*][*][*]:numeric, \  ; Vertical wind   [Pa/s]
                                        opt:logical)
local a, R, cp, k, pi, d2r, name_tim, name_prs, name_lat, name_lon, tim, prs, \
      lat, lon, units, unit, dTdt, dTdx, dTdy, dTdp, dx, dy, logp, Tadv_x,    \
      Tadv_y, Tadv, Tadv_p, Q1
begin
    ; Set constant variable
    a   = 6.37122e06        ; Radius of the earth                 [m]
    R   = 287.              ; Individual gas constant for dry air [J/kg/K]
    cp  = 1005.             ; Specific heat at constant pressure  [J/kg/K]
    k   = R/cp
    pi  = 4*atan(1.)        ; = 3.141592
    d2r = pi/180.           ; Degree to radian
    
    ; Get coordinate variable
    name_tim  = T!0
    name_prs  = T!1
    name_lat  = T!2
    name_lon  = T!3
    tim       = T&$name_tim$
    prs       = T&$name_prs$
    lat       = T&$name_lat$
    lon       = T&$name_lon$

    ; Change time unit (given -> [s])
    units = tim@units
    unit  = stringtochar(units)
    unit := str_upper(chartostring(unit(0:2)))
    ; if (unit .eq. "mon") then dt = ??? end if    
    if (unit .eq. "DAY") then dt = 86400. end if    
    if (unit .eq. "HOU") then dt =  3600. end if    
    if (unit .eq. "MIN") then dt =    60. end if    
    if (unit .eq. "SEC") then dt =     1. end if    
    print("Warning: Forward/Backward scheme is used for time boundary.")

    ; Change pressure unit (hPa -> Pa)
    if isatt(opt, "prs_unit") then 
        if (opt@prs_unit .eq. "hPa") then 
            prs := prs*100.     ; [Pa]
            print("Note! Converting pressure unit : [hPa]->[Pa]")
        else if (opt@prs_unit .eq. "Pa") then 
            print("Note! Using pressure unit : [Pa]")
        end if
        end if
    else
        print("Fatal: You must set opt@prs_unit.")
        exit
        ; print("Note! If pressure unit is hPa, define opt@prs_unit=hPa.")
    end if

    ; Calculate apparent heat source equation
    dTdt = center_finite_diff_n(T, tim*dt, False, 1, 0)
    dTdx = dTdt
    do iy = 0, dimsizes(lat)-1
        dx = lon*d2r*a*cos(lat(iy)*d2r)
        dTdx(:,:,iy,:) = center_finite_diff_n(T(:,:,iy,:), dx, True, 1, 2)
    end do
    dy   = lat*d2r*a
    dTdy = center_finite_diff_n(T, dy, False, 1, 2)

    ; Finite difference in log-P (d/dp -> 1/p*d/dlogp.)
    logp = log(prs)
    dTdp = center_finite_diff_n(T, logp, False, 1, 1)
    dTdp = dTdp/conform(dTdp, prs, 1)          ; take it back to pressure

    ; Horizontal advection term
    Tadv_x = U*dTdx
    Tadv_y = V*dTdy
    Tadv = Tadv_x + Tadv_y

    ; Vertical advection term (THETA)
    Tadv_p = -W*(k*T/conform(T, prs, 1) - dTdp)

    Q1 = cp*(dTdt + Tadv + Tadv_p)  ; Apparent Heat Source
    copy_VarCoords(T, Q1)
    Q1@long_name = "Apparent Heat Source"
    Q1@units = "J/kg/s"
    
    return Q1
end ; /*}}}*/

; Calculate apparent moisture sink (Q2) /*{{{*/
undef("apparent_moisture_sink") 
function apparent_moisture_sink(q[*][*][*][*]:numeric, \  ; Specific humidity [kg/kg]
                                U[*][*][*][*]:numeric, \  ; Zonal wind        [m/s]  
                                V[*][*][*][*]:numeric, \  ; Meridional wind   [m/s]
                                W[*][*][*][*]:numeric, \  ; Vertical wind     [Pa/s] 
                                          opt:logical)
local a, L, pi, d2r, name_tim, name_prs, name_lat, name_lon, tim, prs, lat, lon, \
      units, unit, dqdt, dqdx, dqdy, dqdp, dx, dy, logp, qadv_x, qadv_y, qadv,   \
      qadv_p, Q2
begin
    ; Set constant variable
    a   = 6.37122e06        ; Radius of the earth                  [m]
    L   = 2.501e6           ; Latent heat of vaporization at 273 K [J/kg]
    pi  = 4*atan(1.)        ; = 3.141592
    d2r = pi/180.           ; Degree to radian

    ; Get coordinate variable
    name_tim  = q!0
    name_prs  = q!1
    name_lat  = q!2
    name_lon  = q!3
    tim       = q&$name_tim$
    prs       = q&$name_prs$
    lat       = q&$name_lat$
    lon       = q&$name_lon$

    ; Change time unit (given -> [s])
    units = tim@units
    unit  = stringtochar(units)
    unit := str_upper(chartostring(unit(0:2)))
    ; if (unit .eq. "mon") then dt = ??? end if    
    if (unit .eq. "DAY") then dt = 86400. end if    
    if (unit .eq. "HOU") then dt =  3600. end if    
    if (unit .eq. "MIN") then dt =    60. end if    
    if (unit .eq. "SEC") then dt =     1. end if    
    print("Warning: Forward/Backward scheme is used for time boundary.")

    ; Change pressure unit (hPa -> Pa)
    if isatt(opt, "prs_unit") then 
        if (opt@prs_unit .eq. "hPa") then 
            prs := prs*100.     ; [Pa]
            print("Note! Converting pressure unit : [hPa]->[Pa]")
        else if (opt@prs_unit .eq. "Pa") then 
            print("Note! Using pressure unit : [Pa]")
        end if
        end if
    else
        print("Fatal: You must set opt@prs_unit.")
        exit
        ; print("Note! If pressure unit is hPa, define opt@prs_unit=hPa.")
    end if

    ; Calculate apparent moisture sink equation
    dqdt = center_finite_diff_n(q, tim*dt, False, 1, 0)
    dqdx = dqdt
    do iy = 0, dimsizes(lat)-1
        dx = lon*d2r*a*cos(lat(iy)*d2r)
        dqdx(:,:,iy,:) = center_finite_diff_n(q(:,:,iy,:), dx, True, 1, 2)
    end do
    dy   = lat*d2r*a
    dqdy = center_finite_diff_n(q, dy, False, 1, 2)

    ; Finite difference in log-P (d/dp -> 1/p*d/dlogp.)
    logp = log(prs)
    dqdp = center_finite_diff_n(q, logp, False, 1, 1)
    dqdp = dqdp/conform(dqdp, prs, 1)          ; take it back to pressure

    ; Horizontal advection term
    qadv_x = U*dqdx
    qadv_y = V*dqdy
    qadv   = qadv_x + qadv_y
    
    ; Vertical advection term
    qadv_p = W*dqdp

    Q2 = -L*(dqdt + qadv + qadv_p)  ; Apparent Moisture Sink
    copy_VarCoords(q, Q2)
    Q2@long_name = "Apparent Moisture Sink"
    Q2@units = "J/kg/s"
    
    return Q2
end ; /*}}}*/

