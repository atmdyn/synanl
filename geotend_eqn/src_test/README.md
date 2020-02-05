## Code explanation
* `poisson2d.ncl` : Poisson 2d Solver
* `poisson3d.ncl` : Poisson 3d Solver
* `plot_poisson3d.ncl` : Ploting Poisson 3d using `poisson3d.ncl` result
* `verify_poisson3d.ncl` : Verify (calculate forcing)
* `shaeC.ncl` : Test code for spherical harmonic synthesis 
* `laplacian.ncl` : Test code (laplacian operator)
* `test_laplace.ncl` : Test code (laplacian operator)
  
You can test above scripts in Rossby server.  
```
$ ncl poisson3d.ncl          # make NetCDF file
$ ncl plot_poisson3d.ncl     # Draw result
$ ncl verify_poisson3d.ncl   # for verify
```  

## How to solve Poisson 3d equation?
X, Y 축에 대해서는 구면 조화 방법을, Z 축에 대해서는 SOR 방법을 이용해 계산한다. 여기서 풀고자하는 3차원 푸아송 방정식은 다음과 같다.  
<!-- 3D Poisson's equation -->
<!-- \frac { dT }{ dt } ={ \triangledown  }^{ 2 }+F -->
![equation](https://latex.codecogs.com/gif.latex?\frac&space;{&space;dT&space;}{&space;dt&space;}&space;={&space;\triangledown&space;}^{&space;2&space;}&plus;F)  
여기서 F는 강제력을, T는 온도를 의미한다. 위의 수식에서 차분법을 적용하면

<!-- 시간 차분을 적용한 방정식 -->
<!-- \frac { { T }_{ i+1 }-{ T }_{ i } }{ \Delta t } ={ \triangledown  }^{ 2 }+\quad F -->
![equation](https://latex.codecogs.com/gif.latex?\frac&space;{&space;{&space;T&space;}_{&space;i&plus;1&space;}-{&space;T&space;}_{&space;i&space;}&space;}{&space;\Delta&space;t&space;}&space;={&space;\triangledown&space;}^{&space;2&space;}&plus;\quad&space;F)  

<!-- { T }_{ i+1 }\quad =\quad { T }_{ i }\quad +\quad \Delta t\quad \left( { \nabla  }^{ 2 }{ T }_{ i }\quad +\quad F \right) -->
![equation](https://latex.codecogs.com/gif.latex?{&space;T&space;}_{&space;i&plus;1&space;}\quad&space;=\quad&space;{&space;T&space;}_{&space;i&space;}\quad&space;&plus;\quad&space;\Delta&space;t\quad&space;\left(&space;{&space;\nabla&space;}^{&space;2&space;}{&space;T&space;}_{&space;i&space;}\quad&space;&plus;\quad&space;F&space;\right))
의 수식을 얻을 수 있다. 여기서 라플라스 연산자에 대한 항에 대해 구면조화를 적용하면(Ref. equation 15 참고),  

<!-- 라플라스 연산자 항에 대한 결과 -->
<!-- ({ \triangledown  }^{ 2 }{ T }_{ i }+F)\quad =\quad ({ { \triangledown  }^{ 2 } }+\frac { { { f }_{ 0 } }^{ 2 } }{ { N }^{ 2 } } \frac { { \partial  }^{ 2 } }{ { \partial  }z^{ 2 } } )\sum { \sum { { a }_{ lm }(z){ Y }_{ lm } }  } +\quad \sum { \sum { b_{ lm }(z){ Y }_{ lm } }  } -->
![equation](https://latex.codecogs.com/gif.latex?({&space;\triangledown&space;}^{&space;2&space;}{&space;T&space;}_{&space;i&space;}&plus;F)\quad&space;=\quad&space;({&space;{&space;\triangledown&space;}^{&space;2&space;}&space;}&plus;\frac&space;{&space;{&space;{&space;f&space;}_{&space;0&space;}&space;}^{&space;2&space;}&space;}{&space;{&space;N&space;}^{&space;2&space;}&space;}&space;\frac&space;{&space;{&space;\partial&space;}^{&space;2&space;}&space;}{&space;{&space;\partial&space;}z^{&space;2&space;}&space;}&space;)\sum&space;{&space;\sum&space;{&space;{&space;a&space;}_{&space;lm&space;}(z){&space;Y&space;}_{&space;lm&space;}&space;}&space;}&space;&plus;\quad&space;\sum&space;{&space;\sum&space;{&space;b_{&space;lm&space;}(z){&space;Y&space;}_{&space;lm&space;}&space;}&space;})  
의 수식을 얻으며 이를 정리하면  

<!-- \sum { \sum { [\quad -\frac { l(l+1) }{ { r }^{ 2 } } { a }_{ lm }(z){ \quad +\quad b }_{ lm }(z)\quad +\quad  }  } \frac { { { f }_{ 0 } }^{ 2 } }{ { N }^{ 2 } } \frac { { a }_{ lm }(z+1)+{ a }_{ lm }(z-1)-2{ a }_{ lm }(z) }{ \Delta { Z }^{ 2 } } ]{ Y }_{ lm } -->
![equation](https://latex.codecogs.com/gif.latex?\sum&space;{&space;\sum&space;{&space;[\quad&space;-\frac&space;{&space;l(l&plus;1)&space;}{&space;{&space;r&space;}^{&space;2&space;}&space;}&space;{&space;a&space;}_{&space;lm&space;}(z){&space;\quad&space;&plus;\quad&space;b&space;}_{&space;lm&space;}(z)\quad&space;&plus;\quad&space;}&space;}&space;\frac&space;{&space;{&space;{&space;f&space;}_{&space;0&space;}&space;}^{&space;2&space;}&space;}{&space;{&space;N&space;}^{&space;2&space;}&space;}&space;\frac&space;{&space;{&space;a&space;}_{&space;lm&space;}(z&plus;1)&plus;{&space;a&space;}_{&space;lm&space;}(z-1)-2{&space;a&space;}_{&space;lm&space;}(z)&space;}{&space;\Delta&space;{&space;Z&space;}^{&space;2&space;}&space;}&space;]{&space;Y&space;}_{&space;lm&space;})  
의 결과를 얻는다. 이를 푸아송 방정식에 대입한 뒤, 정리하면

<!-- X, Y에 대해 구면조화를 적용한 푸아송 방정식 -->
<!-- \sum { \sum { \left( { a }_{ lm }^{ i+1 }\quad =\quad { a }_{ lm }^{ i }\left( z \right) \quad -\quad \frac { { f }_{ 0 }^{ 2 } }{ { N }^{ 2 } } \frac { 2\Delta t }{ { \Delta z }^{ 2 } } { a }_{ lm }^{ i }\left( z \right) \quad -\quad \frac { l\left( l+1 \right)  }{ r } \Delta t\quad { a }_{ lm }^{ i }\left( z \right) \quad +\quad \frac { { f }_{ 0 }^{ 2 } }{ { N }^{ 2 } } \frac { \Delta t }{ { \Delta z }^{ 2 } } \left( { a }_{ lm }^{ i }\left( z+1 \right) +{ a }_{ lm }^{ i }\left( z-1 \right)  \right) \quad +\quad \Delta t\quad { b }_{ lm }^{ i }\left( z \right)  \right) { Y }_{ lm } }  } -->
![equation](https://latex.codecogs.com/gif.latex?\sum&space;{&space;\sum&space;{&space;\left(&space;{&space;a&space;}_{&space;lm&space;}^{&space;i&plus;1&space;}\quad&space;=\quad&space;{&space;a&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z&space;\right)&space;\quad&space;-\quad&space;\frac&space;{&space;{&space;f&space;}_{&space;0&space;}^{&space;2&space;}&space;}{&space;{&space;N&space;}^{&space;2&space;}&space;}&space;\frac&space;{&space;2\Delta&space;t&space;}{&space;{&space;\Delta&space;z&space;}^{&space;2&space;}&space;}&space;{&space;a&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z&space;\right)&space;\quad&space;-\quad&space;\frac&space;{&space;l\left(&space;l&plus;1&space;\right)&space;}{&space;r&space;}&space;\Delta&space;t\quad&space;{&space;a&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z&space;\right)&space;\quad&space;&plus;\quad&space;\frac&space;{&space;{&space;f&space;}_{&space;0&space;}^{&space;2&space;}&space;}{&space;{&space;N&space;}^{&space;2&space;}&space;}&space;\frac&space;{&space;\Delta&space;t&space;}{&space;{&space;\Delta&space;z&space;}^{&space;2&space;}&space;}&space;\left(&space;{&space;a&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z&plus;1&space;\right)&space;&plus;{&space;a&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z-1&space;\right)&space;\right)&space;\quad&space;\\&space;&plus;\quad&space;\Delta&space;t\quad&space;{&space;b&space;}_{&space;lm&space;}^{&space;i&space;}\left(&space;z&space;\right)&space;\right)&space;{&space;Y&space;}_{&space;lm&space;}&space;}&space;})  
를 얻는다. 이는 X, Y 축에 대해 구면 조화를 적용한 3차원 푸아송 방정식이며 Z축에 대해서는 구면 조화를 적용하기 까다로우므로, 여기서는 반복법을 이용해 푼다. 여기서는  

<!-- relaxation (w part) -->
<!-- relaxation\quad =\quad \frac { { { f }_{ 0 } }^{ 2 } }{ { N }^{ 2 } } \frac { 2\Delta t }{ \Delta { Z }^{ 2 } } -->
![equation](https://latex.codecogs.com/gif.latex?relaxation\quad&space;=\quad&space;\frac&space;{&space;{&space;{&space;f&space;}_{&space;0&space;}&space;}^{&space;2&space;}&space;}{&space;{&space;N&space;}^{&space;2&space;}&space;}&space;\frac&space;{&space;2\Delta&space;t&space;}{&space;\Delta&space;{&space;Z&space;}^{&space;2&space;}&space;})  
의 값을 조정하여 해에 수렴하는 속도를 조절할 수 있다.  


## Reference
* spherical harmonics - [link](http://scipp.ucsc.edu/~haber/ph116C/SphericalHarmonics_12.pdf)

