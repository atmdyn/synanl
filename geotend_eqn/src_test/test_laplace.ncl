n = 42  ; T42
l = 1
Ts = new((/2,n,n/), "float")
Ts = 0.
Ts(0,l,1) = 1.  ; l1 m1
nlon = 180

set_sphere_radius(1.)
R = get_sphere_radius()
print("R = "+R)

; Sphrical harmonics synthesis
Te = shseC(Ts, nlon)     ; Fixed grid
; Te = shsgC(Ts, nlon)   ; Gaussian grid

lapT = lapsF(Te)
printVarSummary(lapT)

Ts2 = shaeC(lapT)
printVarSummary(Ts2)

print(Ts2(0,l,1))               ; return -l*(l+1)
print(1./sqrt(-Ts2(0,l,1)/(l*(l+1))))  ; return R
