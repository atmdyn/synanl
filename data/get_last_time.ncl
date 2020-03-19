; Getting time index after 1000 days.

; Read data
begidx = 4000
f = addfile("atmos_daily.nc", "r")
ps = f->ps(begidx:,:,:)
ucomp = f->ucomp(begidx:,:,:,:)
vcomp = f->vcomp(begidx:,:,:,:)
temp  = f->temp(begidx:,:,:,:)
omega = f->omega(begidx:,:,:,:)
vor   = f->vor(begidx:,:,:,:)
height    = f->height(begidx:,:,:,:)
tdt_ndamp = f->tdt_ndamp(begidx:,:,:,:)
lonb  = f->lonb
latb  = f->latb
phalf = f->phalf

; Output data
setfileoption("nc","Format","LargeFile")
fc = addfile("atmos_daily_1000-last.nc", "c")
fc->ps = ps
fc->ucomp     = ucomp
fc->vcomp     = vcomp
fc->temp      = temp 
fc->omega     = omega
fc->vor       = vor  
fc->height    = height
fc->tdt_ndamp = tdt_ndamp
fc->lonb  = lonb 
fc->latb  = latb 
fc->phalf = phalf
