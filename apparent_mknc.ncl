; ------------------------------------------------ ;
; Created on  Thu Jan 02  18:37:01 2020
; @author : Hyeong-Gyu Kim.
;
;-- Description
;  - Calculate Q1(apparent heat source) and Q2(apparent moisture sink)
;   using function in "./apparent_calc.ncl"
;  - Input variable must be 4 dimension([Time]x[plev]x[lat]x[lon])
;  - In this code, [Pa] of pressure unit is used.
;  - Note! You must set opt@prs_unit
;    ex) opt@prs_unit = "hPa"
;        opt@prs_unit = "Pa"
;  - You can run this code in Rossby server.
; 
;-- Update
;  - 2020-01-02 : Maked by ADLab
;
;-- Reference
;  - Yanai et al., (1973)
;  - Lee et al., (2008)
;  - Fueglistaler et al., (2009)
; ------------------------------------------------ ;

load "./apparent_calc.ncl"

y4m2 = "201607"
path = "/data2/eraint/1.5x1.5/"

; Read file
fU = addfile(path+"Upl/"+"Upl."+y4m2+".6hr.grib", "r")
U  = fU->U_GDS0_ISBL 
fV = addfile(path+"Vpl/"+"Vpl."+y4m2+".6hr.grib", "r")
V  = fV->V_GDS0_ISBL 
fW = addfile(path+"Wpl/"+"Wpl."+y4m2+".6hr.grib", "r")
W  = fW->W_GDS0_ISBL 
fT = addfile(path+"Tpl/"+"Tpl."+y4m2+".6hr.grib", "r")
T  = fT->T_GDS0_ISBL 
fQ = addfile(path+"Qpl/"+"Qpl."+y4m2+".6hr.grib", "r")
q  = fQ->Q_GDS0_ISBL

opt = True
opt@prs_unit = "hPa"    ; Pressure unit in input data is "hPa"
Q1 = apparent_heat_source(T, U, V, W, opt)
Q2 = apparent_moisture_sink(q, U, V, W, opt)

setfileoption("nc","Format","LargeFile")
ofile = "apparent_heat_moist.nc"
system("rm -f "+ofile)
fc = addfile(ofile, "c")
fc->Q1=Q1
fc->Q2=Q2
