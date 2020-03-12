load "/home/gue/work/TODO/synopanal/apparent_eqn/apparent_calc.ncl"

y4m2 = "201601"
path = "/data2/eraint/1.5x1.5/"
fin  = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
U    = fin->U_GDS0_ISBL(:,{200:1000},:,:)
fin  = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
V    = fin->V_GDS0_ISBL(:,{200:1000},:,:)
fin  = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T    = fin->T_GDS0_ISBL(:,{200:1000},:,:)    ; Temperature [K]
fin  = addfile(path+"Wpl/"+"Wpl."+y4m2+".6hr.grib", "r")
W    = fin->W_GDS0_ISBL(:,{200:1000},:,:)
fQ   = addfile(path+"Qpl/"+"Qpl."+y4m2+".6hr.grib", "r")
q    = fQ->Q_GDS0_ISBL(:,{200:1000},:,:)

opt = True
opt@prs_unit = "hPa"

Q1 = apparent_heat_source(T,U,V,W,opt)
Q2 = apparent_moisture_sink(q,U,V,W,opt)

setfileoption("nc","Format","LargeFile")
ofile = "apparent_heat_moist.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->Q1=Q1
fc->Q2=Q2
