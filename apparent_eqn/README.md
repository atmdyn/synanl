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
![Apparent_heat_source](http://www.sciweavers.org/tex2img.php?eq=%5Cfrac%20%7B%20%5Cpartial%20T%20%7D%7B%20%5Cpartial%20t%20%7D%20%5Cquad%20%2B%5Cquad%20v%5Ccdot%20gradT%5Cquad%20-%5Cquad%20%5Comega%20%5Cleft%28%20%5Cfrac%20%7B%20%5Ckappa%20T%20%7D%7B%20p%20%7D%20-%5Cfrac%20%7B%20%5Cpartial%20T%20%7D%7B%20%5Cpartial%20p%20%7D%20%20%5Cright%29%20%5Cquad%20%3D%5Cquad%20%5Cfrac%20%7B%20Q%20%7D%7B%20%7B%20c%20%7D_%7B%20p%20%7D%20%7D%20&bc=White&fc=Black&im=jpg&fs=12&ff=arev&edit=0)

* Lee et al., (2008) - JMSJ
![Apparent_moisgure_sink]()


(Q1: Apparent heat source [J/kg/s], Q2: Apparent moisture sink [K/s])

<p align="center"><img width="50%" src="/images/apparent.png" /></p>
![apparent_reuslt](/images/apparent.png)   <!-- refer to Emanuel (1994) pp 501 -->

  
## Reference
* Yanai, M., Esbensen, S., & Chu, J. H. (1973). Determination of bulk properties of tropical cloud clusters from large-scale heat and moisture budgets. Journal of the Atmospheric Sciences, 30(4), 611-627. [link](https://journals.ametsoc.org/doi/abs/10.1175/1520-0469(1973)030%3C0611:DOBPOT%3E2.0.CO;2) 

* Fueglistaler, S., B. Legras, A. Beljaars, J. J. Morcrette, A. Simmons, A. M. Tompkins, and S. Uppala, 2009: The diabatic heat budget of the upper troposphere and lower/mid stratosphere in ECMWF reanalyses. Q. J. R. Meteorol. Soc., 135, 21–37 [link](https://rmets.onlinelibrary.wiley.com/doi/abs/10.1002/qj.361)  

* Lee, D.-K., J.-G. Park, and J. Kim, 2009: Heavy Rainfall Events Lasting 18 Days from July 31 to August 17, 1998, over Korea. J. Meteorol. Soc. Japan, 86, 313–333. [link](https://www.jstage.jst.go.jp/article/jmsj/86/2/86_2_313/_article/-char/ja/)

