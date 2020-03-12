f  = new((/181, 360/), "float")
lat = fspan(-90, 90, 181)
pi  = 4*atan(1.)
d2r = pi/180.
OMG = 7.2921e-5
f  = conform(f, 2*OMG*sin(lat*d2r), 0)
sin2 = sin(lat*d2r)*sin(lat*d2r)
sin2(90) = 1.
f2 = f*f

vs = shaeC(f2)
; vs = tri_trunC(vs,42)   ; T42
; v2 = shseC(vs,240)

; print(vs(0,:,0)+"  "+vs(1,:,0))
vs = 0.

Y00 = (1./3.)*sqrt(4.*pi)
Y20 = (2./3.)*sqrt(4.*pi/5.)
; vs(0,0,0) = Y00
vs(0,0,0) = 1
; vs(0,2,0) = Y20
v2  = shseC(vs,360)
dum = 1./sqrt(4*pi)
print(v2(:,0)/dum)
