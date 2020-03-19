; Read data
fin1 = addfile("atmos_average.nc", "r")
pk   = fin1->pk(:)         ; [phalf]
bk   = fin1->bk(:)         ; [phalf]
zsfc = fin1->zsurf(:,:)    ; [lat] x [lon]

fin2 = addfile("atmos_daily_1000-last.nc", "r")
ps     = fin2->ps       (:,  :,:) ; [time] x [lat] x [lon]
ucomp  = fin2->ucomp    (:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
vcomp  = fin2->vcomp    (:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
temp   = fin2->temp     (:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
omega  = fin2->omega    (:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
height = fin2->height   (:,:,:,:) ; [time] x [pfull] x [lat] x [lon]
tdt    = fin2->tdt_ndamp(:,:,:,:) ; [time] x [pfull] x [lat] x [lon]

dims = dimsizes(ucomp)

p_half = new((/dims(0),dims(1)+1,dims(2),dims(3)/), "float")
p_full = new(dims, "float")
copy_VarCoords(temp(:,0,:,:), p_half(:,0,:,:))
copy_VarCoords(pk, p_half(0,:,0,0))
copy_VarCoords(temp, p_full)

; calculate p_full
do k = 0, dims(1)
    ; p = pk + bk*p_surf (See spectral_dynamics.html)
    p_half(:,k,:,:) = tofloat(pk(k) + bk(k)*ps(:,:,:))/100.
end do

do k = 0, dims(1)-1
    alpha = 1.0 - p_half(:,k,:,:) \
                * ( log(p_half(:,k+1,:,:)) - log(p_half(:,k,:,:)) ) \
                / (     p_half(:,k+1,:,:)  -     p_half(:,k,:,:)  )
    p_full(:,k,:,:) = exp( log(p_half(:,k+1,:,:)) - alpha )
end do    
; For bk, pk = 0, See press_and_geopot.f90, 156 lines
p_full(:,0,:,:) = exp( log(p_half(:,1,:,:)) -1 ) 

; vertical interpolation
pout   = ispan(50,1000,50)
linlog = -2  ; log, extrapolation
U   = int2p_n_Wrap(p_full, ucomp , pout, linlog, 1)
V   = int2p_n_Wrap(p_full, vcomp , pout, linlog, 1)
T   = int2p_n_Wrap(p_full, temp  , pout, linlog, 1)
W   = int2p_n_Wrap(p_full, omega , pout, linlog, 1)
Z   = int2p_n_Wrap(p_full, height, pout, linlog, 1)
linlog = 2  ; log, Note! heating -> no extrapolation
Tdt = int2p_n_Wrap(p_full, tdt   , pout, linlog, 1)

fc1 = addfile("fms_HS_U.nc", "c")
fc1->U = U

fc2 = addfile("fms_HS_V.nc", "c")
fc2->V = V

fc3 = addfile("fms_HS_T.nc", "c")
fc3->T = T

fc4 = addfile("fms_HS_W.nc", "c")
fc4->W = W

fc5 = addfile("fms_HS_Z.nc", "c")
fc5->Z = Z

fc6 = addfile("fms_HS_Tdt.nc", "c")
fc6->Tdt = Tdt

