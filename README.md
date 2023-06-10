# cubic_transformation_rs
 a method to generate a univariate non-normal random variable given first four moments

# Reference

* A.I. Fleishman, “A method for simulating nonnormal distributions,” Psychometrika, vol. 43, pp. 521–532, 1978.
* Kjetil Hoyland, Michal Kaut, and Stein W. Wallace, "A Heuristic for Moment-Matching Scenario Generation," Computational Optimization and Applications, Vol. 24, pp. 169-185, 2003.

# 實作細節

1. 給定一個隨機變數Y的平均值(mean)、變異數(variance)、偏度(skewness)、超峰度(ex-kurtosis)，由公式(不偏或是漸近不偏)反推出前4階(原始)動差(moments)EY1, EY2, EY3, EY4 。
2. 由標準常態分佈取出樣本n_scenario個隨機變數X，且算出X的前12階原始動差 EX1, EX2,..., EX12。
3. cubic transform: Y=a + bX + cX**2 + dX**3, 求出a,b,c,d即可由X生成目標樣本Y。

