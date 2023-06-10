# cubic_transformation_rs
 a method to generate a univariate non-normal random variable given first four moments

# Reference

* A.I. Fleishman, “A method for simulating nonnormal distributions,” Psychometrika, vol. 43, pp. 521–532, 1978.
* Kjetil Hoyland, Michal Kaut, and Stein W. Wallace, "A Heuristic for Moment-Matching Scenario Generation," Computational Optimization and Applications, Vol. 24, pp. 169-185, 2003.

# crate

* [gomez](https://docs.rs/gomez/latest/gomez/index.html)
* [levenberg_marquardt](https://docs.rs/levenberg-marquardt/latest/levenberg_marquardt/)

# 演算法

由標準常態分佈X~N(0,1) 抽出S個樣本。

可由cubic transform Y=a+ bX + cX^2 + dX^3得到所要的隨機樣本Y，
其中Y的mean, variance, skew, kurtosis為指定之值。因此問題在於如何得到參數a,b,c,d之值。

可用nonlinear least-square算法求出a,b,c,d之值, 但必須給定X的前12階動差, error function, 
error function的gradient vector與hessian matrix。

注意a,b,c,d不一定有解，此時給出最小化錯誤平方之參數值即可。

已知關係式 Y=a+ bX + cX^2 + dX^3
兩側取期望值得 E(Y) = a + bE(X) + cE(x^2) + dE(X^3) = f(a,b,c,d)

Y^2 = (a+bX+cX^2 + d X^3)^2  兩側取期望值得E(Y^2)
Y^3 = (a+bX+cX^2 + d X^3)^3  兩側取期望值得E(Y^3)
Y^4 = (a+bX+cX^2 + d X^3)^4  兩側取期望值得E(Y^4)

由於已經給定目標前4階(中央)動差，即可反推出E(Y), E(Y^2), E(Y^3), E(Y^4)
可得least square的目標函數為最小化g=
    (E(Y)- f(a,b,c,d))^2 +
    (E(Y^2)- f^2(a,b,c,d))^2 +
    (E(Y^3)- f^3(a,b,c,d))^2 +
    (E(Y^4)- f^4(a,b,c,d))^2的(a,b,c,d)之值。


可得gradient vector為[


]
 

