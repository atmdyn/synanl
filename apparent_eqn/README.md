## Code explanation
We're going to write code using NCL. Each code is described below
* `apparent_calc.ncl`       : Main code for calculate equation about apparent.
* `apparent_mknc.ncl`       : Make NetCDF file using function in `apparent_calc.ncl`
* `apparent_plot.ncl`       : Draw Q1, Q2 for verify.
* `apparent_script_org.ncl` : Exercise script. (Not essential)  

You can test above scripts in Rossby server.
```
$ ncl apparent_mknc.ncl     # make NetCDF file
$ ncl apparent_plot.ncl     # Draw result
```  
  
## Theoretical framework
* Fueglistaler et al., (2009) - QJR  
![Apparent_heat_source](/images/apparent_eq1.png)  
> <p align="center"><img width="75%" src="/images/apparent.png" /></p>
  
* Lee et al., (2008) - JMSJ  
![Apparent_moisgure_sink](/images/apparent_eq2.png)  
(Q1: Apparent heat source [J/kg/s], Q2: Apparent moisture sink [K/s])

<!-- refer to Emanuel (1994) pp 501 -->
<p align="center"><img width="75%" src="/images/apparent.png" /></p>
<center>Fig. 15.3 Averge values of Q1, Q2, adiabatic cooling, and radiative cooling during phase 3 of GATE.</center>

  
## Reference
* Yanai, M., Esbensen, S., & Chu, J. H. (1973). Determination of bulk properties of tropical cloud clusters from large-scale heat and moisture budgets. Journal of the Atmospheric Sciences, 30(4), 611-627. [link](https://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1973)030%3C0611:DOBPOT%3E2.0.CO;2) 

* Fueglistaler, S., B. Legras, A. Beljaars, J. J. Morcrette, A. Simmons, A. M. Tompkins, and S. Uppala, 2009: The diabatic heat budget of the upper troposphere and lower/mid stratosphere in ECMWF reanalyses. Q. J. R. Meteorol. Soc., 135, 21–37 [link](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.361)  

* Lee, D.-K., J.-G. Park, and J. Kim, 2009: Heavy Rainfall Events Lasting 18 Days from July 31 to August 17, 1998, over Korea. J. Meteorol. Soc. Japan, 86, 313–333. [link](https://www.jstage.jst.go.jp/article/jmsj/86/2/86_2_313/_article/-char/ja/)

