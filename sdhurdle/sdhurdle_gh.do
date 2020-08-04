**************************************
* Stochastic dynamic hurdle program  *
* - Unobserved heterogentity: 		 *
*    Gauss hermit quadrature		 *
**************************************
* Kunz & Winkelmann. 2016. 			 *
* An Econometric Model of Healthcare * 
* Demand With Nonlinear Pricing. 	 *
* Health Economics. 				 *
**************************************


/* Program sth */

/* Program sth */
cap program drop sdhgh
program sdhgh, eclass sortpreserve
	syntax varlist [if] [in] [ , QPoints(integer 12) *]
	version 13.1
	gettoken lhs rhs : varlist
	global qp=`qpoints'
	* starting values
	tempname b0 b1
	qui cloglog `lhs' `rhs'
	mat `b0'= e(b)
	qui ztpnm `lhs' `rhs' if `lhs'>0
	mat `b1'= e(b)
	* regression
	ml model lf0 sdhgh_lf (`lhs'=`rhs') (`rhs') /lnalpha , technique(bfgs)  title(Stochastic dhurdle gauss-hermit regression)
	ml init `b0' `b1'  , copy
	ml max, diff 
end


cap program drop sdhgh_lf
program sdhgh_lf
	version 13.1
	args todo b lnfj
	tempvar lin01 lin11 lin02 lin12 lin03 lin13 lin04 lin14
	tempvar l0 l1 lin0 lin1 y prj k s pr delta spr lnsig sig d0 d1  
	tempvar f1 pd1 y1 lj S ss_pd
	tempvar f2 pd2 y2	
	qui{
	mleval `lin0'  = `b' , eq(1)
	mleval `lin1'  = `b' , eq(2)
	mleval `lnsig' = `b' , eq(3) sca
	local y2 $ML_y1
	g double  `l0' = .
	g double  `l1' = .	
	g double `pd2' = 0
	g double  `f2' = 0
	g double `spr' = 0
	g double `ss_pd' = 0
	
	/* Get points and weights for Gaussian-Hermite quadrature. */
	tempvar x w 
	g double `x'=.
	g double `w'=.
	_GetQuad, avar(`x') wvar(`w') quad($qp)
	replace `x'=`x'*sqrt(2)
	local r=1
	while `r' <=$qp {
			sca `sig' = exp(`lnsig')*`x'[`r']
			replace `l0'    = exp(`lin0'+`sig')
			replace `l1'    = exp(`lin1'+`sig')
			su `y2'
			replace 	`f2' = gammap(`y2',(`l1'-`l0'))
			replace     `f2' = 1                         															if `l1'<`l0'
			forvalues      j = 0/`r(max)' {
				replace `f2' = `f2'-exp(-(`l1'-`l0'))*((`l1'-`l0'))^`j'/lngamma(`j' + 1)  							if `j'<=`y2'-1 & `l1'<`l0'
				}
			replace    `pd2' = poissonp(`l1', `y2')		    														if `l1'==`l0'					//poisson ll
			replace    `pd2' = exp(-`l0')			   																if `l1'!=`l0' & `y2' == 0		//poisson ll (0)
			replace    `pd2' = `l0'*`l1'^(`y2'-1)*exp(-(`l0'))*(`f2'/((`l1'-`l0')^`y2'))							if `l1'!=`l0' & `y2' > 0		//Dhurdle part
		replace `ss_pd'=`ss_pd'+`pd2'*`w'[`r']
		local r =`r'+1
		}
	replace `lnfj'  =  ln(`ss_pd'/sqrt(_pi))
	}
end

