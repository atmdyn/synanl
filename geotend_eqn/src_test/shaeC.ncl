a=addfile("/data2/eraint/1.5x1.5/VOpl/VOpl.201201.6hr.grib","r")

vo = a->VO_GDS0_ISBL(:,{900},:,:)
vs = shaeC(vo)
vs = tri_trunC(vs,42)   ; T42
v2 = shseC(vs,240)
copy_VarCoords(vo,v2)


s = "pdf"
s@xmin = 0
s@xmax = 360
s@ymin = -90
s@ymax =  90
w = mapbase(s,"shaeC_tri_test")
r = True
r@style = "fill"
r@cmax = 2.e-4
r@cmin =-2.e-4
r@cint = 1.e-5
addcontour(w,vo(0,:,:),r)
frame(w)
addcontour(w,v2(0,:,:),r)
