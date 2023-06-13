# cubic_transformation_rs
 a method to generate a univariate non-normal random variable given first four moments

# Reference

* A.I. Fleishman, “A method for simulating nonnormal distributions,” Psychometrika, vol. 43, pp. 521–532, 1978.
* Kjetil Hoyland, Michal Kaut, and Stein W. Wallace, "A Heuristic for Moment-Matching Scenario Generation," Computational Optimization and Applications, Vol. 24, pp. 169-185, 2003.

# crate

* [gomez](https://docs.rs/gomez/latest/gomez/index.html)
* [levenberg_marquardt](https://docs.rs/levenberg-marquardt/latest/levenberg_marquardt/)

# 演算法

給定目標分佈的平均值(mean), 變異數(variance)，偏度(skewness)，超峰度(ex-kurtosis)。
首先由標準常態分佈X~N(0,1) 抽出S個樣本。且計算出樣本X的前12階(原始)動差(moments)。

設定轉換後的樣本隨機變數Z的統計量為[0, 1, skewness，ex-kurtosis+3].
使用least square error algorithm得到Z=a+ bX + cX^2 + dX^3的4個參數(a,b,c,d)。
使用上述公式由隨機樣本X得到隨機樣本Z，此時Z的統計量應接近於[0, 1, skewness，ex-kurtosis+3]。

再使用Y=mean + variance.sqrt() * Z得到最終樣本Y, 此時樣本Y的(mean, variance, skewness, ex-kurtosis)會近似於目標分佈的統計量。
