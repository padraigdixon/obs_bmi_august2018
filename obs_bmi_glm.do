



cd "...data\"

set showbaselevels all //this makes the base category explicit where used.

use data.dta, clear

//Keep only completed episodes

tab epistat

distinct n_eid

/////////////////////////////////////////////////////////////////////////////

*Creating some new variables to support analysis


///


*Waist-hip ratio

gen whr=waist_circum/hip_circum





gen obese_whr=0
replace obese_whr=1 if whr>.85 & sex==0 //female
replace obese_whr=1 if whr>.9 & sex==1 //male  

//Some exclusions

drop if bmi_trad==. & bmi_imp==.


//Consider relationship between bmi-traditional and bmi-impedance

concord bmi_trad bmi_imp if cost_per!=., level(.99)

/*


drop if difference is more than 5 st dev away from mean difference

*/

gen bmi_diff=bmi_trad-bmi_imp

qui su bmi_diff
gen mean_bmi_diff=r(mean)
gen mean_bmi_diff_sd=r(sd)

drop if bmi_diff < mean_bmi_diff-5*mean_bmi_diff_sd
distinct n_eid
drop if bmi_diff>mean_bmi_diff+5*mean_bmi_diff_sd & bmi_diff!=.
distinct n_eid


*Generate a compsite weight variable for analysis, using impedance where traditional measure is absent

gen bmi_combined=.
replace bmi_combined =bmi_trad
replace bmi_combined=bmi_imp if bmi_trad==.



egen bmi_who_comb=cut(bmi_comb), at(10 18.5 20(2.5)30 35 40 80) label //would be good to label with more informative names in due course




gen obese_combined=0
replace obese_comb=1 if bmi_comb>=30


*Generate quintiles of deprivation
xtile quin_deprivation=townsend,n(5) 

*For sensitivity analysis...generate a variable for individuals with no pre-existing conditions. Coding is 1 if *no* pre_existing conditions

gen no_existing_condition=0
replace no_exist=1 if number_self_conditions==0 & number_self_cancers==0

/////

*Create a binary cost variable to use when comparing private health insurance

 gen nonzero_cost=1 if cost_per>0 & cost_per<.
 recode nonzero_cost (.=0)
 
 recode private (-3=.) (-1=.)
 
 tab private
 
 logit nonzero_cost i.private if private!=. ,or
 margins private, level(99)
 marginsplot
 
////////////////////////////////////////////////////////////////////////////
*Descriptive statistics

////////////////////////////////////////////////////////////////////////////



global covariates i.sex i.days_exercise i.alcohol_freq i.qualif i.employment i.quin age_at_recruitment

foreach var of varlist days_exercise alcohol_freq qualif employment quin {

replace `var'=. if `var'==-3|`var'==-1 //Recode prefer not to answer and don't know to missing

}

//On a variable-by-variable basis, recoded "none of the above" to a separate "other" category:
*Qualifications

replace qualifications =8 if qualifications==-7 

*Employment

replace employment=8 if employment==-7



//Covariate missingness

foreach var of varlist sex days_exercise alcohol_freq health_rating qualif employment quin age_at_recruitment {

distinct n_eid if `var'==.

cap gen `var'_miss=0
replace  `var'_miss=1 if `var'==.

su `var'_miss

}

//Summary table of missingness

su *_miss




drop *_miss 



save bmi_obs_analysis200818, replace

use bmi_obs_analysis200818, clear



su bmi_comb whr bn.i.sex bn.i.days_exercise bn.i.alcohol_freq bn.i.qualif bn.i.employment bn.i.quin age_at_recruitment if cost_per!=., sep(0)

bysort bmi_who_comb: su bmi_comb whr bn.i.sex bn.i.days_exercise bn.i.alcohol_freq bn.i.qualif bn.i.employment bn.i.quin age_at_recruitment if cost_per!=., sep(0)






//Exploration of the distribution of costs


kdensity cost_per if cost_per>0


*Density of OLS residuals 
qui reg cost_per bmi_comb $covariates

predict exphat
centile exphat, centile(2 5 6 7)
predict u, residual

kdensity u




//Tests of family and link functions 



foreach var of varlist bmi_comb bmi_who obese_comb {



xi: boxcox cost_per `var' $covariates if cost_per>0 

}



foreach var of varlist bmi_comb bmi_who obese_comb {


qui glm cost_per `var' $covariates,  link(log) family(gamma)
linktest, link(log) family(gamma)


}


foreach var of varlist bmi_comb bmi_who obese_comb {
*** Run GLM to generate residuals for Park test for family
quietly glm cost_per `var' $covariates, link(log) family(gamma)

*** Generate ln(raw residuals squared) and xbetahat for Park test
predict double rawresid`var', response
generate lnrawresid2`var' = ln(rawresid`var'^2)
predict double xbetahat`var', xb

*** Modified Park test 
regress lnrawresid2`var' xbetahat`var', robust


*

}




////////////////////////////////////////////////////////////////////////////////

*ANALYSIS MODELS

////////////////////////////////////////////////////////////////////////////////
cd "...\data\"

use bmi_obs_analysis200818, clear

//Need to set base levels - set base to most frequent value (25-27.5)

fvset base freq bmi_who


*set an excel sheet to capture results and graphs

putexcel set initial_margins,  replace



*BMI models




foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(base_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix base_results_`var'=r(table)
putexcel Z25=matrix(base_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix base_mchange_`var'=r(table)
putexcel T6=matrix(base_mchange_`var'), names
margins, dydx(`var') level(99) post
est store Base_Model_`var'
matrix margin_`var'=e(b)
scalar define margin_mean=margin_`var'[1,1]
marginsplot
graph export base_`var'.png, replace
 putexcel A25 = picture(base_`var') 
 putexcel A6 = etable

}





foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
quietly sum `var' if e(sample)
local sd = r(sd)
margins, at(`var'=gen(`var' + `sd') ) 
matrix predicted_`var'=r(table)
putexcel M6=matrix(predicted_`var'), names
}




foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=gen(`var') ) 
matrix predicted_`var'=r(table)
putexcel P6=matrix(predicted_`var'), names

}



//Comparison of glm with a crude linear model



foreach var of varlist bmi_comb  {
reg cost_per `var' $covariates,   robust   level (99)
margins, dydx(`var') level(99) 
marginsplot
mchange `var',amount (sd) stat(ci) level(99) 

}

//Categorical models

foreach var of varlist bmi_who  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix base_results_`var'=r(table)
putexcel Z25=matrix(base_results_`var'), names
margins, dydx(`var') level(99) post
est store Base_Model_`var'
marginsplot
graph export base_`var'.png, replace
 putexcel A25 = picture(base_`var') 
 putexcel A6 = etable
}


foreach var of varlist bmi_who  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
margins `var', level (99)
matrix predicted_`var'=r(table)
putexcel P6=matrix(predicted_`var'), names

}


foreach var of varlist  obese_comb  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix base_results_`var'=r(table)
putexcel Z25=matrix(base_results_`var'), names
margins, dydx(`var') level(99) post
est store Base_Model_`var'
marginsplot
graph export base_`var'.png, replace
 putexcel A25 = picture(base_`var') 
 putexcel A6 = etable
}




foreach var of varlist obese_comb  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(0 1))
matrix predicted_`var'=r(table)
putexcel P6=matrix(predicted_`var'), names

}
////////////////////////////////////////////////////////////////////////////////


*WHR models



foreach var of varlist whr  {
putexcel set initial_margins,  sheet(base_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix base_results_`var'=r(table)
putexcel Z25=matrix(base_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix base_mchange_`var'=r(table)
putexcel T6=matrix(base_mchange_`var'), names
mchange `var',delta(0.01) stat(ci) level(99) 
matrix base_mchange_`var'=r(table)
putexcel T16=matrix(base_mchange_`var'), names

margins, dydx(`var') level(99) post
est store Base_Model_`var'
marginsplot
graph export base_`var'.png, replace
 putexcel A25 = picture(base_`var') 
 putexcel A6 = etable
mchange `var',amount (sd) stat(ci) level(99) 
mchange `var',delta(0.01) stat(ci) level(99) 

}



foreach var of varlist whr  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
quietly sum `var'  if e(sample)
local sd = r(sd)
margins, at(`var' =gen(`var'  + `sd') ) 
matrix sd_predicted_`var'=r(table)
putexcel M6=matrix(sd_predicted_`var'), names
}



foreach var of varlist whr  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
local increment=0.01
margins, at(`var' =gen(`var'  + `increment') ) 
matrix inc_predicted_`var'=r(table)
putexcel P6=matrix(inc_predicted_`var'), names
}



foreach var of varlist  obese_whr  {
putexcel set initial_margins,  sheet(base_`var') modify

glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix base_results_`var'=r(table)
putexcel Z25=matrix(base_results_`var'), names
margins, dydx(`var') level(99) post
est store Base_Model_`var'
marginsplot
graph export base_`var'.png, replace
 putexcel A25 = picture(base_`var') 
 putexcel A6 = etable
}


foreach var of varlist obese_whr  {
putexcel set initial_margins,  sheet(base_`var') modify
glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(0 1) ) 
matrix predicted_`var'=r(table)
putexcel P6=matrix(predicted_`var'), names

}




////////////////////////////////////////////////////////////////////////////////
*1. SENSTIVITY ANALYSIS - ADJUSTED BMI FOR WHR AND VICE VERSA


*WHR models



foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(whradj_`var') modify

glm cost_per `var' $covariates whr,  link(log) family(gamma) robust eform level (99)
matrix whradj_results_`var'=r(table)
putexcel Z25=matrix(whradj_results_`var'), names
margins, dydx(`var') level(99) post
est store whradj_Model_`var'
marginsplot
graph export whradj_`var'.png, replace
 putexcel A25 = picture(whradj_`var') 
 putexcel A6 = etable
}

foreach var of varlist bmi_who  {
putexcel set initial_margins,  sheet(whradj_`var') modify

glm cost_per i.`var' $covariates whr,  link(log) family(gamma) robust eform level (99)
matrix whradj_results_`var'=r(table)
putexcel Z25=matrix(whradj_results_`var'), names
margins, dydx(`var') level(99) post
est store whradj_Model_`var'
marginsplot
graph export whradj_`var'.png, replace
 putexcel A25 = picture(whradj_`var') 
 putexcel A6 = etable
}


foreach var of varlist whr  {
putexcel set initial_margins,  sheet(whradj_`var') modify

glm cost_per `var' $covariates bmi_comb,  link(log) family(gamma) robust eform level (99)
matrix bmiadj_results_`var'=r(table)
putexcel Z25=matrix(bmiadj_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix bmiadj_mchange_`var'=r(table)
putexcel T6=matrix(bmiadj_mchange_`var'), names

margins, dydx(`var') level(99) post
est store bmiadj_Model_`var'
marginsplot
graph export bmiadj_`var'.png, replace
 putexcel A25 = picture(bmiadj_`var') 
 putexcel A6 = etable



glm cost_per `var' $covariates bmi_who,  link(log) family(gamma) robust eform level (99)
matrix bmiadj_who_results_`var'=r(table)
putexcel Z25=matrix(bmiadj_who_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix bmiadj_who_mchange_`var'=r(table)
putexcel T16=matrix(bmiadj_who_mchange_`var'), names
margins, dydx(`var') level(99) post
est store bmiadj_who_Model_`var'
marginsplot
graph export bmiadj_who_`var'.png, replace
 putexcel A125 = picture(bmiadj_who_`var') 
 putexcel A16 = etable


}

////////////////////////////////////////////////////////////////////////////////
*2. SENSTIVITY ANALYSIS - DROPPING COVARIATES WITH HIGHER MISSINGNESS

*days exercise (5.4%) and qualifications (1.9%) have more than 1% missingness

global shorter_covariates i.sex  i.alcohol_freq  i.employment i.quin age_at_recruitment


*BMI models

foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(shorter_`var') modify

glm cost_per `var' $shorter_covariates,  link(log) family(gamma) robust eform level (99)
matrix short_results_`var'=r(table)
putexcel Z25=matrix(short_results_`var'), names

mchange `var',amount (sd) stat(ci) level(99) 
matrix short_mchange_`var'=r(table)
putexcel T6=matrix(short_mchange_`var'), names

margins, dydx(`var') level(99) post
est store less_covar_`var'
marginsplot
graph export shorter_`var'.png, replace
 putexcel A25 = picture(shorter_`var') 
 putexcel A6 = etable

}



foreach var of varlist bmi_who obese_comb  {
putexcel set initial_margins,  sheet(shorter_`var') modify

glm cost_per i.`var' $shorter_covariates,  link(log) family(gamma) robust eform level (99)
matrix short_results_`var'=r(table)
putexcel Z25=matrix(short_results_`var'), names
margins, dydx(`var') level(99) post
est store less_covar_`var'
marginsplot
graph export shorter_`var'.png, replace
 putexcel A25 = picture(shorter_`var') 
 putexcel A6 = etable
}

////////////////////////////////////////////////////////////////////////////////


*WHR models



foreach var of varlist whr  {
putexcel set initial_margins,  sheet(shorter_`var') modify

glm cost_per `var' $shorter_covariates,  link(log) family(gamma) robust eform level (99)

matrix short_results_`var'=r(table)
putexcel Z25=matrix(short_results_`var'), names
margins, dydx(`var') level(99) post
est store less_covar_`var'
marginsplot
graph export shorter_`var'.png, replace
 putexcel A25 = picture(shorter_`var') 
 putexcel A6 = etable
mchange `var',amount (sd) stat(ci) level(99) 

}



foreach var of varlist  obese_whr  {
putexcel set initial_margins,  sheet(shorter_`var') modify

glm cost_per i.`var' $shorter_covariates,  link(log) family(gamma) robust eform level (99)
matrix short_results_`var'=r(table)
putexcel Z25=matrix(short_results_`var'), names
margins, dydx(`var') level(99) post
est store less_covar_`var'
marginsplot
graph export shorter_`var'.png, replace
 putexcel A25 = picture(shorter_`var') 
 putexcel A6 = etable
}





////////////////////////////////////////////////////////////////////////////////
*Sensitivity analysis
*3. Sensitivity analysis using two-part model 

//Two-part model, estimated as a sensitivity analysis

*BMI models

foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(tpm_`var') modify

tpm cost_per `var' $covariates,  first(logit) s(glm, fam(gamma) link(log)) robust eform level (99)
matrix tpm_results_`var'=r(table)
putexcel Z25=matrix(tpm_results_`var'), names
margins, dydx(`var') level(99) post
est store tpm_`var'
marginsplot
graph export tpm_`var'.png, replace
 putexcel A25 = picture(tpm_`var') 
 putexcel A6 = etable


}



foreach var of varlist bmi_who obese_comb  {
putexcel set initial_margins,  sheet(shorter_`var') modify

tpm cost_per i.`var' $covariates,  first(logit) s(glm, fam(gamma) link(log)) robust eform level (99)
matrix tpm_results_`var'=r(table)
putexcel Z25=matrix(tpm_results_`var'), names
margins, dydx(`var') level(99) post
est store tpm_`var'
marginsplot
graph export tpm_`var'.png, replace
 putexcel A25 = picture(tpm_`var') 
 putexcel A6 = etable
}

////////////////////////////////////////////////////////////////////////////////


*WHR models



foreach var of varlist whr  {
putexcel set initial_margins,  sheet(shorter_`var') modify

tpm cost_per `var' $covariates,  first(logit) s(glm, fam(gamma) link(log)) robust eform level (99)
matrix tpm_results_`var'=r(table)
putexcel Z25=matrix(tpm_results_`var'), names
margins, dydx(`var') level(99) post
est store tpm_`var'
marginsplot
graph export tpm_`var'.png, replace
 putexcel A25 = picture(tpm_`var') 
 putexcel A6 = etable
mchange `var',amount (sd) stat(ci) level(99)

}



foreach var of varlist  obese_whr  {
putexcel set initial_margins,  sheet(shorter_`var') modify

tpm cost_per i.`var' $covariates,  first(logit) s(glm, fam(gamma) link(log)) robust eform level (99)
matrix tpm_results_`var'=r(table)
putexcel Z25=matrix(tpm_results_`var'), names
margins, dydx(`var') level(99) post
est store tpm_`var'
marginsplot
graph export tpm_`var'.png, replace
 putexcel A25 = picture(tpm_`var') 
 putexcel A6 = etable
}




////////////////////////////////////////////////////////////////////////////////
*Sensitivity analysis
*4. Only never smokers i.e if smoking_status==0 - make sure the preserve/restore commands are actually implemented

preserve

keep if smoking==0

*BMI models

foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(never_smoke_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_smoke_results_`var'=r(table)
putexcel Z25=matrix(no_smoke_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix nosmoke_mchange_`var'=r(table)
putexcel T6=matrix(nosmoke_mchange_`var'), names

margins, dydx(`var') level(99) post
est store no_smoke_`var'
marginsplot
graph export never_smoke_`var'.png, replace
 putexcel A25 = picture(never_smoke_`var') 
 putexcel A6 = etable


}

restore

preserve

keep if smoking==0

foreach var of varlist bmi_who obese_comb  {
putexcel set initial_margins,  sheet(never_smoke_`var') modify

glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_smoke_results_`var'=r(table)
putexcel Z25=matrix(no_smoke_results_`var'), names
margins, dydx(`var') level(99) post
est store no_smoke_`var'
marginsplot
graph export never_smoke_`var'.png, replace
 putexcel A25 = picture(never_smoke_`var') 
 putexcel A6 = etable

}

restore

////////////////////////////////////////////////////////////////////////////////


*WHR models


preserve

keep if smoking==1

foreach var of varlist whr  {
putexcel set initial_margins,  sheet(never_smoke_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix never_smoke_results_`var'=r(table)
putexcel Z25=matrix(never_smoke_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix never_smoke_mchange_`var'=r(table)
putexcel T6=matrix(never_smoke_mchange_`var'), names

mchange `var',amount(0.1) stat(ci) level(99)
matrix never_smoke_mchange_`var'=r(table)
putexcel T16=matrix(never_smoke_mchange_`var'), names

margins, dydx(`var') level(99) post
est store no_smoke_`var'
marginsplot
graph export never_smoke_`var'.png, replace
 putexcel A25 = picture(never_smoke_`var') 
 putexcel A6 = etable


}

restore

preserve

keep if smoking==1

foreach var of varlist  obese_whr  {
putexcel set initial_margins,  sheet(never_smoke_`var') modify

glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix never_smoke_results_`var'=r(table)
putexcel Z25=matrix(never_smoke_results_`var'), names
margins, dydx(`var') level(99) post
est store no_smoke_`var'
marginsplot
graph export never_smoke_`var'.png, replace
 putexcel A25 = picture(never_smoke_`var') 
 putexcel A6 = etable
mchange `var',amount (sd) stat(ci) level(99) 

}


restore


////////////////////////////////////////////////////////////////////////////////
*Sensitivity analysis
*5. Only those without pre-existing health conditions

preserve 

keep if no_existing_condition==1

*BMI models

foreach var of varlist bmi_comb  {
putexcel set initial_margins,  sheet(no_existing_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_exist_results_`var'=r(table)
putexcel Z25=matrix(no_exist_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99) 
matrix no_existing_mchange_`var'=r(table)
putexcel T6=matrix(no_existing_mchange_`var'), names

margins, dydx(`var') level(99) post
est store no_existing_`var'
marginsplot
graph export no_existing_`var'.png, replace
 putexcel A25 = picture(no_existing_`var') 
 putexcel A6 = etable




}

restore

preserve

keep if no_existing_condition==1

foreach var of varlist bmi_who obese_comb  {
putexcel set initial_margins,  sheet(no_existing_`var') modify

glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_exist_results_`var'=r(table)
putexcel Z25=matrix(no_exist_results_`var'), names
margins, dydx(`var') level(99) post
est store no_existing_`var'
marginsplot
graph export no_existing_`var'.png, replace
 putexcel A25 = picture(no_existing_`var') 
 putexcel A6 = etable
}

restore

////////////////////////////////////////////////////////////////////////////////


*WHR models

preserve 

keep if no_existing_condition==1

foreach var of varlist whr  {
putexcel set initial_margins,  sheet(no_existing_`var') modify

glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_exist_results_`var'=r(table)
putexcel Z25=matrix(no_exist_results_`var'), names
mchange `var',amount (sd) stat(ci) level(99)
matrix no_existing_mchange_`var'=r(table)
putexcel T6=matrix(no_existing_mchange_`var'), names

mchange `var',amount(0.1) stat(ci) level(99)
matrix no_existing_mchange_`var'=r(table)
putexcel T16=matrix(no_existing_mchange_`var'), names


margins, dydx(`var') level(99) post
est store no_existing_`var'
marginsplot
graph export no_existing_`var'.png, replace
 putexcel A25 = picture(no_existing_`var') 
 putexcel A6 = etable


}

restore

preserve 

keep if no_existing_condition==1

foreach var of varlist  obese_whr  {
putexcel set initial_margins,  sheet(no_existing_`var') modify

glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
matrix no_exist_results_`var'=r(table)
putexcel Z25=matrix(no_exist_results_`var'), names
margins, dydx(`var') level(99) post
est store no_existing_`var'
marginsplot
graph export no_existing_`var'.png, replace
 putexcel A25 = picture(no_existing_`var') 
 putexcel A6 = etable
}


restore

//Close the excel file used to capture results

putexcel close

////////////////////////////////////////////////////////////////////////////////

*Some graphs for the paper

*1. Continuous BMI at representative values

foreach var of varlist bmi_comb  {
qui glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(bmi=(12(1)70)) level(99) 
marginsplot, saving(cont_bmi, replace) recast(line)  ciopts(recast(rline)lpattern(dash)) 

}



*2. Categorical BMI


foreach var of varlist bmi_who  {
qui glm cost_per i.`var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, dydx(`var') level(99) 
marginsplot, saving(cat_bmi) recast(line)  ciopts(recast(rline)lpattern(dash))

}

*3. Age, stratified by sex

foreach var of varlist bmi_comb  {
qui glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(age=(40(5)75) sex=(0 1)) level(99) 
marginsplot, saving(age_sex, replace) recast(line)  ciopts(recast(rline)lpattern(dash))


}



////////////////////////////////////////////////////////////////////////////////



//Coefplot to summarise sensitivity analyses
  


foreach var of varlist bmi_comb  {
 glm cost_per `var' $covariates,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store Base_`var'
}


foreach var of varlist bmi_comb  {
glm cost_per `var' $covariates whr,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store WHR_`var'
}



foreach var of varlist bmi_comb  {
glm cost_per `var' $shorter_covariates if smoking==0,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store Nosmoke_`var'

}

foreach var of varlist bmi_comb  {
glm cost_per `var' $shorter_covariates,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store Shortcovar_`var'

}

foreach var of varlist bmi_comb  {

twopm cost_per `var' $covariates,  first(logit) s(glm, fam(gamma) link(log)) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store twopm_`var'

}




foreach var of varlist  bmi_comb  {
glm cost_per `var' $covariates if no_existing_condition==1,  link(log) family(gamma) robust eform level (99)
margins, at(`var'=(12(1)70)) level(99) post
est store No_exist_`var'
}



//Code to make continuous graph 

coefplot (Base_bmi_combined, label("Base model")) (Test_2_bmi_combined,label("WHR-adjusted BMI")) (Nosmoke_bmi_combined,label("Never smokers")) (Shortcovar_bmi_combined,label("Limited covariates")) /// 
(No_exist_bmi_combined,label("No existing conditions")) (twopm_bmi_combined,label("2-part model")), at ytitle("Predicted Cost") xtitle("BMI")   /// 
recast(line)lwidth(*2)ciopts(recast(rline)lpattern(dash)) saving(sensitivity_cont_bmi)



//Code to make graph of marginal effects 

coefplot (Base_Model_bmi_combined \ whradj_Model_bmi_combined \  no_smoke_bmi_combined \  less_covar_bmi_combined  \ no_existing_bmi_combined \ tpm_bmi_combined),  mlabels aseq swapnames xline(`=margin_mean')  saving("marginal effect by model")


