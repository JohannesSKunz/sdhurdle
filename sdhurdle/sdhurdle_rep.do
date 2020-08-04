**************************************
* Stochastic dynamic hurdle program  *
* - Reporting time missmatch		 *
**************************************
* Kunz & Winkelmann. 2016. 			 *
* An Econometric Model of Healthcare * 
* Demand With Nonlinear Pricing. 	 *
* Health Economics. 				 *
**************************************

/* Program dhurdlerep */
cap program drop sdhrep
program sdhrep, eclass sortpreserve
	syntax varlist [if] [in]
	version 13.1
	gettoken lhs rhs : varlist
	// starting values
	tempname b0 b1 b2 
	qui sdh `lhs' `rhs'
	mat `b2'= e(b)
	// Dhurdle regression
	ml model lf0 sdhrep_lf (`lhs'=`rhs') (`rhs') , tech(nr)  title(Stochastic dhurdel reporting time regression)
	ml init `b2'  , copy
	ml max , diff 
end


/* Program sth_lf */
cap program drop sdhrep_lf
program sdhrep_lf
	version 13.1
	args todo b lnfj
	tempvar lin01 lin11 lin02 lin12 lin03 lin13 lin04 lin14
	tempvar l0 l1 lin0 lin1 y prj k s pr delta spr lnsig sig d0 d1  
	tempvar f1 pd1 y1 lj S ss_pd
	tempvar f2 pd2 y2	
		qui{
		mleval `lin0'  = `b' , eq(1)
		mleval `lin1'  = `b' , eq(2)
		local y  $ML_y1
		g double  `l0' = .
		g double  `l1' = .	
		g double `prj' = 0
		g double `pd1' = 0
		g double `pd2' = 0
		g double  `y1' = 0
		g double  `y2' = 0
		g double  `f1' = 0
		g double  `f2' = 0
		g double  `pr' = 0
		g double `spr' = 0
		g double `ss_pd' = 0
		g double  `d0' = delta
		g double  `d1' = 1-`d0'
		
			replace `l0'    = exp(`lin0')
			replace `l1'    = exp(`lin1')
			su `y'
			local k=r(max)
			forvalues  s  = 0/`k' {
				//first prob
				replace     `y1' = `s'  
				su `y1'
				replace     `f1' = gammap(`y1',(`l1'-`l0')*`d1')
				replace     `f1' = 1                         																if `l1'<`l0'
				forvalues      j = 0/`r(max)' {
					replace `f1' = `f1'-exp(-(`l1'-`l0')*`d1')*((`l1'-`l0')*`d1')^`j'/round(exp(lnfactorial(`j')),1)  		if `j'<=`y1'-1 & `l1'<`l0'
					}
				replace    `pd1' = poissonp(`l1'*`d1', `y1')		    										    		if `l1'==`l0'					//poisson ll
				replace    `pd1' = `l0'*`l1'^(`y1'-1)*exp(-(`l0'*`d1'))*(`f1'/((`l1'-`l0')^`y1'))							if `l1'!=`l0' & `y1' > 0		//Dhurdle part
				replace    `pd1' = exp(-`l0'*`d0')*`pd1'		  +(1-exp(-`l0'*`d0'))*poissonp(`l1'*`d1',`y1')   			if `l1'!=`l0' 
				replace    `pd1' = exp(-`l0'*`d0')*exp(-`l0'*`d1')+(1-exp(-`l0'*`d0'))*exp(-`l1'*`d1')						if `l1'!=`l0' & `y1' == 0		//poisson ll (0)
				//second prob
				replace     `y2' = `y'-`s' 
				su `y2'
				replace 	`f2' = gammap(`y2',(`l1'-`l0')*`d0')
				replace     `f2' = 1                         															   	if `l1'<`l0'
				forvalues      j = 0/`r(max)' {
					replace `f2' = `f2'-exp(-(`l1'-`l0')*`d0')*((`l1'-`l0')*`d0')^`j'/round(exp(lnfactorial(`j')),1)  		if `j'<=`y2'-1 & `l1'<`l0'
					}
				replace    `pd2' = poissonp(`l1'*`d0', `y2')		    													if `l1'==`l0'					//poisson ll
				replace    `pd2' = exp(-`l0'*`d0')			   																if `l1'!=`l0' & `y2' == 0		//poisson ll (0)
				replace    `pd2' = `l0'*`l1'^(`y2'-1)*exp(-(`l0'*`d0'))*(`f2'/((`l1'-`l0')^`y2'))							if `l1'!=`l0' & `y2' > 0			//Dhurdle part
				//sum up 
				replace    `prj' = `prj'+ `pd1'*`pd2'   if `pd2'!=. & `pd1'!=. & `s'<=`y'
				}
	replace `lnfj'  =  ln(`prj')
	}
end
