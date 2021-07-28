********************************************************************************
* COMP_STAT  Create comparative-statics figures.
********************************************************************************


*** load data
* loop through sheets: BETA, RHO1, RHO2, ALPHA, GAMMA, NU
foreach sheet in BETA RHO1 RHO2 ALPHA GAMMA NU {
	
	* import
	import excel using "${DIR_SIMULATIONS}/sim_comp_stat.xlsx", sheet("`sheet'") firstrow case(lower) clear
	
	* rename and label key variables
	local var_name = strlower("`sheet'")
	rename par_value `var_name'
	if "`sheet'" == "BETA" label var beta "Discount factor, {&beta}"
	else if "`sheet'" == "RHO1" label var rho1 "Infection parameter, {&rho}{subscript:1}"
	else if "`sheet'" == "RHO2" label var rho2 "Infection parameter, {&rho}{subscript:2}"
	else if "`sheet'" == "ALPHA" label var alpha "Cobb-Douglas share of intermediate inputs, {&alpha}"
	else if "`sheet'" == "GAMMA" label var gamma "Relative productivity of infected workers, {&gamma}"
	else if "`sheet'" == "NU" label var nu "Value of life, {&nu}"
	
	* save
	save "${DIR_SIMULATIONS}/`var_name'.dta", replace
}

* clear memory
clear

* loop through files: BETA, RHO1, RHO2, ALPHA, GAMMA, NU
foreach sheet in BETA RHO1 RHO2 ALPHA GAMMA NU {
	
	* append
	local var_name = strlower("`sheet'")
	append using "${DIR_SIMULATIONS}/`var_name'.dta"
}


*** format data
* rename weekly discount factor
rename beta beta_weekly
label var beta_weekly "Weekly discount factor, {&beta}"

* create annual discount factor
gen double beta_annual = beta_weekly^52
label var beta_annual "Annual discount factor, {&beta}{superscript:52}"

* normalize GDP to = 1 with no pandemic
foreach var in gdp_commit gdp_no_commit gdp_no_lock gdp_no_pand {
	replace `var' = `var'/gdp_no_pand
}

* convert death rate into percent
foreach var in deaths_commit deaths_no_commit deaths_no_lock deaths_no_pand {
	replace `var' = `var'*100
}

* normalize value of planner's problem to = 1 with no pandemic
foreach var in v_commit v_no_commit v_no_lock v_no_pand {
	replace `var' = `var'/v_no_pand
}

* save
order beta_weekly beta_annual rho1 rho2 alpha gamma nu gdp_no_pand gdp_commit gdp_no_commit gdp_no_lock deaths_no_pand deaths_commit deaths_no_commit deaths_no_lock v_no_pand v_commit v_no_commit v_no_lock is_baseline
save "${DIR_SIMULATIONS}/sim_comp_stat.dta", replace



*** plot figures
* extract baseline x-axis values
foreach var of varlist beta_annual rho1 rho2 alpha gamma nu {
	sum `var' if `var' < . & is_baseline == 1, meanonly
	global `var'_baseline = r(mean)
}

* plot GDP and death rate vs. annual discount factor, beta
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit beta_annual if beta_annual < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${beta_annual_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_beta, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_beta.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit beta_annual if beta_annual < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${beta_annual_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_beta, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_beta.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit beta_annual if beta_annual < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${beta_annual_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.99 (.01)1.02, grid gstyle(dot) gmin gmax format(%4.3f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_beta, replace)
if $fig_export graph export "${DIR_FIGURES}/value_beta.${fig_format}", replace

* plot GDP and death rate vs. infection parameter, rho1
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit rho1 if rho1 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho1_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.1f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_rho1, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_rho1.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit rho1 if rho1 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho1_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.1f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_rho1, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_rho1.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit rho1 if rho1 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho1_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.99 (.01)1.02, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_rho1, replace)
if $fig_export graph export "${DIR_FIGURES}/value_rho1.${fig_format}", replace

* plot GDP and death rate vs. infection parameter, rho2
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit rho2 if rho2 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho2_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.1f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_rho2, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_rho2.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit rho2 if rho2 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho2_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.1f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_rho2, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_rho2.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit rho2 if rho2 < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${rho2_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.99 (.01)1.02, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_rho2, replace)
if $fig_export graph export "${DIR_FIGURES}/value_rho2.${fig_format}", replace

* plot GDP and death rate vs. Cobb-Douglas share of intermediate inputs, alpha
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit alpha if alpha < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${alpha_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%2.1f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_alpha, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_alpha.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit alpha if alpha < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${alpha_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%2.1f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_alpha, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_alpha.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit alpha if alpha < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${alpha_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.9(.4)1.7, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_alpha, replace)
if $fig_export graph export "${DIR_FIGURES}/value_alpha.${fig_format}", replace

* plot GDP and death rate vs. relative productivity of infected workers, gamma
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit gamma if gamma < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${gamma_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%2.1f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_gamma, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_gamma.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit gamma if gamma < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${gamma_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%2.1f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_gamma, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_gamma.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit gamma if gamma < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${gamma_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(0.985(.010)1.015, grid gstyle(dot) gmin gmax format(%4.3f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_gamma, replace)
if $fig_export graph export "${DIR_FIGURES}/value_gamma.${fig_format}", replace

* plot GDP and death rate vs. value of life, nu
tw ///
	(connected gdp_no_pand gdp_no_lock gdp_commit gdp_no_commit nu if nu < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${nu_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%1.0f)) ylabel(.7(.1)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("GDP (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(GDP_nu, replace)
if $fig_export graph export "${DIR_FIGURES}/GDP_nu.${fig_format}", replace

tw ///
	(connected deaths_no_pand deaths_no_lock deaths_commit deaths_no_commit nu if nu < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${nu_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%1.0f)) ylabel(0(.7)2.1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	ytitle("Death rate (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(deaths_nu, replace)
if $fig_export graph export "${DIR_FIGURES}/deaths_nu.${fig_format}", replace

tw ///
	(connected v_no_pand v_no_lock v_commit v_no_commit nu if nu < ., sort lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${nu_baseline}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ylabel(.97(.04)1.09, grid gstyle(dot) gmin gmax format(%4.3f)) ///
	ytitle("Value (rel. to no pandemic)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	ysize(2.8) xsize(7) scale(1.9) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) ///
	name(value_nu, replace)
if $fig_export graph export "${DIR_FIGURES}/value_nu.${fig_format}", replace
