load "./poisson_solver.ncl"

fin = addfile("/data2/eraint/1.5x1.5/Tpl/Tpl.201607.6hr.grib", "r")
T0 = fin->T_GDS0_ISBL(2,{200:1000},{35},{127})


; Wavelength
Lz = 200.     ; [hPa]
Ly = 10.
Lx = 10. 

; Coordinate
lev = T0&lv_ISBL1
lat = ispan(-90,  90, 1)    ; 181
lon = ispan(  0, 359, 1)    ; 360
nlev = dimsizes(lev)
nlat = dimsizes(lat)
nlon = dimsizes(lon)

;-- Forcing Part
F = new((/nlev,nlat,nlon/), "float")
F!0  = "lev"
F!1  = "lat"
F!2  = "lon"
F&lev = lev
F&lat = lat
F&lon = lon
hgt3 = conform(F, lev, 0)
lat3 = conform(F, lat, 1)
lon3 = conform(F, lon, 2)

A = 10./86400.   ; amplitude 10 [K/day]
F = A* exp(- ( ( lat3        / Ly )^2 \
             + ((lon3-180)   / Lx )^2 \    ; 90
             + ((hgt3-500) / Lz )^2 ) )
opt = True
opt@T0 = T0
opt@lev_units = "hPa"
opt@dt = 3.e8
opt@thres = 1.e-2
; opt@upper_bound = 
; opt@lower_bound = 

solve = poisson_solver(F, opt)

ofile = "run_output.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->forcing  = F
fc->solution = solve
