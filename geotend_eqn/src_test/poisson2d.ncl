; Poisson's equation solver (test version)

; Coordinate
lat = ispan(-90, 90, 1)
lon = ispan(0, 359, 1)

Lx = 10.
Ly = 10.
; set_sphere_radius(1.)
R = get_sphere_radius()

f1 = new((/181,360/), "float") ; forcing
lat2 = conform(f1, lat, 0)
lon2 = conform(f1, lon, 1)

; Gaussian forcing (90: warming, 270: cooling)
f1 = exp(-((lat2/Ly)^2 + ((lon2-90)/Lx)^2))       ; 90
f1 = f1 - exp(-((lat2/Ly)^2 + ((lon2-270)/Lx)^2)) ; 270

; Te     = f1*0 + 273.15
; lapT   = lapsF(Te)
; Ts_LHS = shaeC(lapT)    ; -(l*(l+1))/R^2 a_lm Y_lm
Ts_RHS = shaeC(-f1)       ; b_lm Y_lm
; printVarSummary(Ts_RHS)   ; [2] x [181] x [181]

a_lm = todouble(Ts_RHS)*0.
do l = 1, dimsizes(lat)-1
    a_lm(:,l,:) = -(R^2)/(l*(l+1))*Ts_RHS(:,l,:)
end do    
; printVarSummary(a_lm)

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
; res@cnLinesOn = False
plt = gsn_csm_contour(wks, f1, res)
plt = gsn_csm_contour(wks, a_lm(0,:,:), res)
plt = gsn_csm_contour(wks, shseC(a_lm, dimsizes(lon)), res)
