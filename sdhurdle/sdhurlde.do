**************************************
* Stochastic dynamic hurdle program  *
**************************************
* Kunz & Winkelmann. 2016. 			 *
* An Econometric Model of Healthcare * 
* Demand With Nonlinear Pricing. 	 *
* Health Economics. 				 *
**************************************

/* Program sth */
cap program drop sdhurdle
program sdh, eclass sortpreserve
	syntax varlist [if] [in]
	gettoken lhs rhs : varlist
	version 13.1
	* starting values
	tempname b0 b1
	qui cloglog `lhs' `rhs'
	mat `b0'= e(b)
	qui ztp `lhs' `rhs' if `lhs'>0
	mat `b1'= e(b)
	* regression
	ml model lf sdh_lf (`lhs'=`rhs') (`rhs'), technique(nr) title(Stochastic dhurdel regression)
	ml init `b0' `b1', copy
	ml max
end

/* Program sth_lf */
cap program drop sdh_lf
program sdh_lf
	version 13.1
	args lnf lin0 lin1
tempvar l0 l1 y
	gen double `l0'=exp(`lin0')
	gen double `l1'=exp(`lin1')
	qui gen double `y'=$ML_y1
tempvar f1
	qui sum `y'
	qui gen double `f1'=gammap(`y',`l1'-`l0')
	qui replace `f1'=1 if `l1'<`l0'
	forvalues j = 0/`r(max)' {
	qui replace `f1'=`f1'-exp(-(`l1'-`l0'))*(`l1'-`l0')^`j'/exp(lnfactorial(`j')) if `j'<=`y'-1 & `l1'<`l0'
	}
	/* Calculating log-likelihood */
	qui replace `lnf'=cond(`l1'==`l0' | `y' == 0, ln(poissonp(`l0', `y')), `lin0'+`lin1'*(`y'-1)-`l0'+ln(abs(`f1'))-ln(abs(`l1'-`l0'))*`y')
end

