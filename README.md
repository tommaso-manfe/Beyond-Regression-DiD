# DID-simulations
Montecarlo simulations of difference-in-difference estimating methods.

In this folder you can find an example of the code that I am currently employing in master's thesis work. 
The aim of the paper is to understand the performance of the traditional difference in difference estimator with respect to new semi-parametric methods proposed by the literature when some of its assumptions are relaxed. Following the debiased machine-learning literature, the thesis aims also at evaluating whether the introduction of machine learning first-stage estimates improves the performance of the estimators. 


I consider four different group, C0, C1, T0, T1, which I compute separately to have the possibility to create some heterogeneity in the characteristics of the groups. 

C0 is the control group (d) at t=0 (in the script, it is when d_i=0 and t_i=0)
C1 is the control group at t=1 (in the script, it is when d_i=0 and t_i=1)
T0 is the treated group at t=0 (in the script, it is when d_i=1 and t_i=0)
T1 is the treated group at t=0 (in the script, it is when d_i=1and t_i=1)

For each group, I simulate a multivariate normal distribution of 5 covariates. Contrarily to independent normal random variable draws, taking the multivariate distribution enables to create correlations between the different covariates. I add an unobservable characteristic that is correlated with the probability of being treated (and so also with other covariates). Therefore, without taking first difference, the time invariant unobservable characteristic would create bias in the estimation of the parameter of interest. 

I ensemble the four group in a unique dataset. I create the potential outcomes Y00, Y01, Y10, Y11 (Ydt), which depend on the covariates, the time effect (for t=1), the belonging to the treatment group (for d=1) and the treatment effect (for t=1 and d=1). In addition, I include a random noise component. Then, I create the observed outcome, that is the only potential outcome that can be observed. 

Finally, I artificially create a set of covariates that are available to the researcher but do not describe the data generating process of the dependent variable. This emulates a real circumstance in which the data generating process is sparse (depends on few covariates) but there are many other variables to be used. Such a high-dimensional setting should favour the use of machine learning techniques such as LASSO. The choice of a LASSO model may not only be driven by to basically selects a subset of informative covariates, but also it may be the only way for the estimation in case the number of covariates (and their interaction for example) exceeds the number of observations. When the estimators have first stages on a subset of the data, as in the triple matching estimator, the risk of not having enough observation is higher. I finally estimate the coefficients with the different techniques.

OLS refers to linear regression model that does not exploit the time-dimension of the data.
DID refers to the standard difference-in-differences model that uses linear regression.
IPWDID refers to inverse probability weighting (IPW) model proposed by Abadie (2005). The code is the one implemented by Sant’Anna and Zhao (2020).
ORDID refers to the outcome regression (OR) model of Heckman, Ichimura and Todd (1997). The code is the one implemented by Sant’Anna and Zhao (2020).
IMPDRDID is the doubly robust model of Sant’Anna and Zhao (2020).
TRIPLEMATCHING is derived by the authors’ independent work.
IMPDRDID is the doubly robust model of Sant’Anna and Zhao (2020).
DID_AMLE is the Augmented Minimax Linear Estimation (AMLE) difference in differences model proposed by Lu, Nie and Wager (2019). The model is doubly robust and allows for machine learning first stages estimation.

REFERENCES
Abadie, A. (2005), “Semiparametric difference-in-difference estimators,” Review of Economic Studies, 72, 1-19.
Heckman, J. J., Ichimura, H., and Todd, P. (1997), “Matching as an econometric evaluation estimator: Evidence from
evaluating a job training programme,” The Review of Economic Studies, 64(4), 605–654.
Sant’Anna and Zhao (2020), 'Doubly robust difference-in-differences estimators' Journal of Econometrics
Volume 219, Issue 1, November 2020, Pages 101-122.
Lu, Nie and Wager (2019), Nonparametric Heterogeneous Treatment Effect Estimation in Repeated Cross Sectional Designs, Arkiv.

