#!/bin/csh

set odir = "/work2/hgkim/2022_gdaps/era5/clim"
set vars = ( "Zpl" "Tpl" "Upl" "Vpl" )
set months = ( "01" "02" "03" "04" "05" "06" "07" "08" "09" "10" "11" "12" )


# Main
foreach var ( $vars )
    if ($var == "Zpl") set gvar  = "var129"
    if ($var == "Zpl") set vname = 'zclm'
    if ($var == "Tpl") set gvar  = "var130"
    if ($var == "Tpl") set vname = 'tclm'
    if ($var == "Upl") set gvar  = "var131"
    if ($var == "Upl") set vname = 'uclm'
    if ($var == "Vpl") set gvar  = "var132"
    if ($var == "Vpl") set vname = 'vclm'

    foreach mon ( $months )
        rm -f ${odir}/lndata/*

        # Link 30 years data for climatology (1982 - 2021)
        ln -s /data3/era5/1.5x1.5/${var}/${var}.198[2-9]${mon}.6hr.grib ${odir}/lndata/.
        ln -s /data3/era5/1.5x1.5/${var}/${var}.199[0-9]${mon}.6hr.grib ${odir}/lndata/.
        ln -s /data3/era5/1.5x1.5/${var}/${var}.201[0-9]${mon}.6hr.grib ${odir}/lndata/.
        ln -s /data3/era5/1.5x1.5/${var}/${var}.202[0-1]${mon}.6hr.grib ${odir}/lndata/.
 
        set flis = `ls -1 ${odir}/lndata/*`
        if ( ${#flis} != 30 ) then
            echo "Some file is missing!"
        endif

        cdo -f nc copy -chname,${gvar},${vname} -yhouravg -mergetime ${flis} ${odir}/${var}/dummy.nc
        if ( ${mon} == "02" ) then  # delete 2/29
            cdo delete,month=2,day=29 ${odir}/${var}/dummy.nc ${odir}/${var}/dummy2.nc
            rm -f ${odir}/${var}/dummy.nc
            mv -f ${odir}/${var}/dummy2.nc ${odir}/${var}/dummy.nc
        endif
        cdo -yhouravg ${odir}/${var}/dummy.nc ${odir}/${var}/${vname}"_1982-2021_"${mon}.nc
        rm -f ${odir}/${var}/dummy.nc
    end
end
