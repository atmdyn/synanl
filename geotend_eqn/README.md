## Code explanation
* Making Tendency Term (dZ/dt) (working on it...)
* Studying SOR method  
  
  
## SOR method
SOR(Successive Over Relaxaion) 방법론은 주어진 방정식을 반복적으로 구하여 해에 근접한 값을 얻어내는 방법이다. 이와 유사한 방법으로는 자코비(Jacobi) 방법과 가우스-자이델(Gauss-Seidel) 방법론이 있으며, 해에 도달하는 속도에 차이가 있을 뿐 결과는 동일하게 얻어진다는 차이가 있다.  

예를 들어 아래와 같은 2차원의 확산 방정식을 가정하여 보자. 
  
![equation](/images/geotend1.png)  

확산 방정식이 평형 상태에 도달하게 되면 확산이 없어지므로 
![equation](/images/geotend2.png)
이 되며, 결과적으로 라플라스 방정식의 해를 구하는 것과 동일하게 된다. 즉, 반복법을 이용해 라플라스 방정식의 해를 얻으면, 평형 상태의 확산 방정식 해를 얻을 수 있다.  

확산 방정식의 해를 구하기위해 중앙차분법을 적용하며, 문제를 단순화하기 위하여 x, y 방향으로의 격자 크기가 h로 동일하다고 가정하면  
  
![equation](/images/geotend3.png)  
![equation](/images/geotend4.png)  
![equation](/images/geotend5.png)  
로 단순화할 수 으며, 오메가의 값이 1일때 가우스-자이델 방법이, 1보다 클 경우에는 SOR 방법이 적용된다. 
  
## Reference
Alejandro Garcia - Numerical Methods for Physics-Addison-Wesley (1999)
