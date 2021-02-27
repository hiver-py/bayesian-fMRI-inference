# bayesian-fMRI-inference
Comparison of Bayesian and frequentist dual regression on resting state networks using three different models.
1. A frequentist Dual Regression,
2. A Bayesian model with heteroscedastic variance and uncorrelated error terms 
3. A Bayesian model with heteroscedastic variance and correlated error terms. 

Non-in formative prior distributions are used for both Bayesian models. 

These 3 models can account for varying amounts of information in the data due to varying complexity of the covariance structure. Difference was observed in the subject specific maps, indicating that covariance should be included when making inference.


![fMRI](https://user-images.githubusercontent.com/73787550/109378289-9fd58280-78d1-11eb-9b9d-cc4a08b5aaaa.png)

