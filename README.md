# FunctionalLasso
R codes for the article Roche, A. (2023) Lasso in infinite dimension: application to variable selection in functional multivariate linear regression, Electronic Journal of Statistics, 17(2), 3357-3405.

The code is published in the spririt of the article https://www.nature.com/articles/467753a. It does not include necessarily good comments and may not be robust to the treatment of another dataset. Please do not hesitate to send any comment or question to angelina.roche@u-paris.fr. 

Please find here a description of the files: 
* Main files: 
  - main_elec.R: code for the treatment of electric consumption data (the data file name is energydata_complete.csv and comes from the website https://archive.ics.uci.edu/dataset/374/appliances+energy+prediction
  - main.R: code for the simulation study conducted in the article.
 
* Files containing the R functions:
  - glmnet_func2.R: algorithm to approach the solution of the Lasso criterion inspired from the glmnet algorithm (Friedman, Hastie and Tibshirani, 2010) and Yang and Zou (2015).
    
        J. Friedman, T. Hastie, and R. Tibshirani. Regularization paths for generalized linear models via coordinate descent. Journal of Statistical Software,
33(1):1–22, 2010.

      Y. Yang and H. Zou. A fast unified algorithm for solving group-lasso penalize learning problems. Stat. Comput., 25(6):1129–1141, 2015.

      
  - normes_ps_L2.R: functions allowing to approach norms and scalar products for multivariate functional data
  - tikhonov.R: stochastic descent gradient algorithm to approach the solution of the Ridge criterion for multivariate functional data.
  - fonctions_simuX.R: functions for the simulation of functional data. 
