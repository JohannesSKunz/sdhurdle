This folder contains Monte Carlo simulations and estimation programs for the
stochastic dynamic hurdle model developed in the paper: 
	Kunz & Winkelmann, 2016, An Econometric Model of Health Care Demand With Nonlinear Pricing, Health Economics, forthcoming. 

Note: the ztpnm package needs to be installed:  
	Farbmacher, H. (2011). “Estimation of hurdle models for overdispersed count data” Stata Journal. 11, Number 1, pp. 82–94.

Version: 01 Oct 2016
__________________________________________________________________________________
Programs contain:

sdhurlde
- Estimates the simple stochastic dynamic hurdle model.

sdhurdle_gh
- Accounts for unobserved heterogeneity Gauss-Hermit-Quadrature with 12 points default, might be changed by: “, qp(30) “.

sdhurdle_rep
-Accounts for reporting time mismatch, delta (distance to beginning of period)
has to be defined in the data set (see example below, more information can be found in the paper).

sdhurdle_rep_gh
- Accounts for both, unobserved heterogeneity and reporting mismatch.

The programs are stored separately and have to be loaded in the beginning of
the program.

__________________________________________________________________________________
Example Data

Small (random) SOEP excerpt, containing number of doctor visits,
gender, log income, and age.

The variable "delta" needs to be provided for the reporting time regressions.
It has been constructed using the syntax, for a non-leap year:
g delta=. ;
replace delta=(daynr)/90        if quarter==1 ;
replace delta=(daynr-90)/91  	if quarter==2 ;
replace delta=(daynr-181)/92   	if quarter==3 ;
replace delta=(daynr-273)/92    if quarter==4 ;

__________________________________________________________________________________
For any questions refer to: johannesskunz@gmail.com
