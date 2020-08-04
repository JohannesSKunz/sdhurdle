clear all 
set more off 
use data.dta
set seed 1
cap log close 
log using example.txt, text replace

//call programm
qui do ../sdhurlde.do
qui do ../sdhurdle_gh.do
qui do ../sdhurdle_rep.do     
qui do ../sdhurdle_rep_gh.do

su 

//estimate
sdh	     nrdoctorvisits age  male linc 
sdhrep   nrdoctorvisits age  male linc 
sdhgh    nrdoctorvisits age  male linc , qp(30)
sdhrepgh nrdoctorvisits age  male linc , qp(30)

cap log close 
