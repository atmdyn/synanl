## Code explanation
We're going to write code using NCL. Each code is described below
* `waf_calc.ncl`       : Main code for calculate wave activity flux (WAF).
* `waf_mknc.ncl`       : Make NetCDF file using function in `waf_calc.ncl`
* `waf_plot.ncl`       : Draw Fx, Fy and geopotential height for verify.
* `run_mknc_clim.csh`  : Make climatology for calculate WAF.

You can test above scripts in Rossby server.
```
$ ./run_mknc_clim.csh       # make climatology using cdo
$ ncl waf_mknc.ncl          # make NetCDF file
$ ncl waf_plot.ncl          # Draw result
```  
  
## Original Program
You can download original script in [here](http://www.atmos.rcast.u-tokyo.ac.jp/nishii/programs/index.html)

  
## Reference
*  Takaya, K., & Nakamura, H. (2001). A formulation of a phase-independent wave-activity flux for stationary and migratory quasigeostrophic eddies on a zonally varying basic flow. Journal of the Atmospheric Sciences, 58(6), 608-627. [link](https://doi.org/10.1175/1520-0469(2001)058%3C0608:AFOAPI%3E2.0.CO;2)

