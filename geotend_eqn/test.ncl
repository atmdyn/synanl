; 3D Poisson's equation solver (test version)

; Coordinate
lev = ispan(0, 20, 1)   ; km
lat = ispan(-90, 90, 1)
lon = ispan(0, 359, 1)
nlon = dimsizes(lon)

Lx = 10.
Ly = 10.
Lz = 2  ; km
; set_sphere_radius(1.)
; toint(), tofloat() make segmentation fault
; R = tofloat(get_sphere_radius()) ;<- error
R = get_sphere_radius()     

f1 = new((/21,181,360/), "float") ; forcing
hgt3 = conform(f1, lev, 0)
lat3 = conform(f1, lat, 1)
lon3 = conform(f1, lon, 2)

; Gaussian forcing (90: warming, 270: cooling)
A = 10/86400.   ; amplitude 10 [K/day]
f1 = exp(-((lat3/Ly)^2 + ((lon3-90)/Lx)^2) + ((hgt3-10)/Lz)^2)  ; 90
f1 = A*( f1 - exp(-((lat3/Ly)^2 + ((lon3-270)/Lx)^2 + ((hgt3-10)/Lz)^2)) ); 270

b_lm = shaeC(f1)
; printVarSummary(b_lm)   ; [2] x [21] x [181] x [181]

; Relaxation coefficient
f0 = 1.e-4
N  = 1.e-2
; relax = (f0/N)^2*(2*dt/dz^2)
relax = 0.1  ; = 1 : normal relaxation
             ; > 1 : over relaxationa
resid = 10. ; critical value
dz = 1000.  ; [m]
dt = ((N/f0)^2*(dz)^2)/2.*relax

a_lm = todouble(b_lm*0.)    ; Because R is double
n  = ispan(0,dimsizes(lat)-1,1)
nn1 = conform(a_lm(:,0,:,:),n*(n+1),1)  ; [2] x [l] x [m]
err = 1.d5

do while (err .gt. resid)
    dummy = a_lm
    do k = 1, dimsizes(lev)-2   ; omit boundary
        laps = (nn1/R^2)*dt*a_lm(:,k,:,:) 
        laps(:,0,0) = 0.        ; No mean value for laplacian
        a_lm(:,k,:,:) = a_lm(:,k,:,:) - relax*a_lm(:,k,:,:) \
                        - laps \
                        + relax/2.*(a_lm(:,k+1,:,:) + a_lm(:,k-1,:,:) ) \
                        + dt*b_lm(:,k,:,:)
    end do    
    err = avg(abs(a_lm - dummy))
    print(err+"")
end do    

solution = shseC(a_lm, nlon)
system("rm -f ./output.nc")
fc = addfile("./output.nc", "c")
fc->solution = solution

wks = gsn_open_wks("x11", "nonamed")
res = True
res@cnFillOn = True
; res@cnLinesOn = False
; plt = gsn_csm_contour(wks, f1, res)
; plt = gsn_csm_contour(wks, a_lm(0,:,:), res)
; plt = gsn_csm_contour(wks, shseC(a_lm, dimsizes(lon)), res)
plt = gsn_csm_contour(wks, solution(0,:,:), res)
