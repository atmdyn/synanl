; make geopotential height using hydro() function

fin1 = addfile("atmos_average.nc", "r")
pk   = fin1->pk(:)         ; [phalf]
bk   = fin1->bk(:)         ; [phalf]
zsfc = fin1->zsurf(:,:)    ; [lat] x [lon]

fin2 = addfile("atmos_daily_1000-last.nc", "r")
; ps   = fin2->ps(:,:,:)     ; [time] x [lat] x [lon]
temp = fin2->temp(:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
pfull = temp&pfull
ps = 100000

p_half = tofloat(pk + bk*ps)/100.
p_full = new((/20/), "float")
do k = 1, 19
    alpha = 1.0 - p_half(k) \
                 * ( log(p_half(k+1)) - log(p_half(k)) ) \
                 / (p_half(k+1) - p_half(k))
    p_full(k) = exp( log(p_half(k+1)) - alpha )
end do  
p_full(0) = exp( log(p_half(1)) -1 ) ; bk, pk = 0, See press_and_geopot.f90, 156 lines

print(pfull+"  "+p_full)
