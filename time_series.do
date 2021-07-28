********************************************************************************
* TIME_SERIES  Create time-series figures.
********************************************************************************


*** load data
* import
import excel using "${DIR_SIMULATIONS}/sim_time_series.xlsx", sheet("Sheet1") firstrow clear


*** format data
* rename and label key variables
rename Time week
label var week "Week"
rename Vaccine vacc
label var vacc "Ind: Vaccine has arrived?"
rename lock0 lockdown_nol
label var lockdown_nol "Lockdown policy with no lockdown"
rename lockC lockdown_commit
label var lockdown_commit "Lockdown policy with commitment"
rename lockN lockdown_nocommit
label var lockdown_nocommit "Lockdown policy with no commitment"
rename Y0 gdp_nol
label var gdp_nol "GDP with no lockdown"
rename YC gdp_commit
label var gdp_commit "GDP with commitment"
rename YN gdp_nocommit
label var gdp_nocommit "GDP with no commitment"
rename X0 inv_nol
label var inv_nol "Investment with no lockdown"
rename XC inv_commit
label var inv_commit "Investment with commitment"
rename XN inv_nocommit
label var inv_nocommit "Investment with no commitment"
rename LOGC_0 util_cons_nol
label var util_cons_nol "Flow utility of consumption with no lockdown"
rename LOGC_C util_cons_commit
label var util_cons_commit "Flow utility of consumption with commitment"
rename LOGC_N util_cons_nocommit
label var util_cons_nocommit "Flow utility of consumption with no commitment"
rename NU_0 util_life_nol
label var util_life_nol "Flow utility of being alive with no lockdown"
rename NU_C util_life_commit
label var util_life_commit "Flow utility of being alive with commitment"
rename NU_N util_life_nocommit
label var util_life_nocommit "Flow utility of being alive with no commitment"
rename S0 s_nol
label var s_nol "Share susceptible (S) with no lockdown"
rename SC s_commit
label var s_commit "Share susceptible (S) with commitment"
rename SN s_nocommit
label var s_nocommit "Share susceptible (S) with no commitment"
rename I0 i_nol
label var i_nol "Share infected (I) with no lockdown"
rename IC i_commit
label var i_commit "Share infected (I) with commitment"
rename IN i_nocommit
label var i_nocommit "Share infected (I) with no commitment"
rename R0 r_nol
label var r_nol "Share recovered (R) with no lockdown"
rename RC r_commit
label var r_commit "Share recovered (R) with commitment"
rename RN r_nocommit
label var r_nocommit "Share recovered (R) with no commitment"
rename D0 d_nol
label var d_nol "Share dead (D) with no lockdown"
rename DC d_commit
label var d_commit "Share dead (D) with commitment"
rename DN d_nocommit
label var d_nocommit "Share dead (D) with no commitment"

* create week = 0 as last week before pandemic begins (i.e., pandemic begins in week 1)
qui count
local N_plus_1 = r(N) + 1
set obs `N_plus_1'
replace week = 0 if week == .

* compute rounded-up maximum week number
global week_plot_max = 70

* compute period when miracle cure arrives
qui count if vacc == 1 & inrange(week, 0, ${week_plot_max})
global period_miracle_cure = ${week_plot_max} - `=r(N)' + 1

* generate share of population that is alive
gen float alive_nol = s_nol + i_nol + r_nol
label var alive_nol "Share alive with no lockdown"
gen float alive_commit = s_commit + i_commit + r_commit
label var alive_commit "Share alive with commitment"
gen float alive_nocommit = s_nocommit + i_nocommit + r_nocommit
label var alive_nocommit "Share alive with no commitment"

* generate population flow utility from consumption
gen float util_cons_pop_nol = alive_nol*util_cons_nol
label var util_cons_pop_nol "Population flow utility from consumption with no pandemic"
gen float util_cons_pop_commit = alive_commit*util_cons_commit
label var util_cons_pop_commit "Population flow utility from consumption with commitment"
gen float util_cons_pop_nocommit = alive_nocommit*util_cons_nocommit
label var util_cons_pop_nocommit "Population flow utility from consumption with no commitment"

* generate population flow utility from being alive
gen float util_life_pop_nol = alive_nol*util_life_nol
label var util_life_pop_nol "Population flow utility from being alive with no pandemic"
gen float util_life_pop_commit = alive_commit*util_life_commit
label var util_life_pop_commit "Population flow utility from being alive with commitment"
gen float util_life_pop_nocommit = alive_nocommit*util_life_nocommit
label var util_life_pop_nocommit "Population flow utility from being alive with no commitment"

* generate total flow utility = utility of consumption + utility of being alive
gen float util_tot_pop_nol = util_cons_pop_nol + util_life_pop_nol
label var util_tot_pop_nol "Total population flow utility with no lockdown"
gen float util_tot_pop_commit = util_cons_pop_commit + util_life_pop_commit
label var util_tot_pop_commit "Total population flow utility with commitment"
gen float util_tot_pop_nocommit = util_cons_pop_nocommit + util_life_pop_nocommit
label var util_tot_pop_nocommit "Total population flow utility with no commitment"

* replace variables in week 0 (i.e., before pandemic begins)
replace vacc = 0 if week == 0
foreach var of varlist lockdown* {
	replace `var' = 0 if week == 0
}
foreach var of varlist s_* alive* {
	replace `var' = 1 if week == 0
}
foreach var of varlist i_* r_* d_* {
	replace `var' = 0 if week == 0
}
sort week
foreach var of varlist gdp* {
	replace `var' = gdp_nol[2] if week == 0
}
foreach var of varlist inv* {
	replace `var' = inv_nol[2] if week == 0
}
foreach var of varlist util_cons_* {
	replace `var' = util_cons_nol[2] if week == 0
}
foreach var of varlist util_cons_pop_* {
	replace `var' = util_cons_pop_nol[2] if week == 0
}
foreach var of varlist util_life_* {
	replace `var' = util_life_pop_nol[2] if week == 0
}
foreach var of varlist util_life_* {
	replace `var' = util_life_pop_nol[2] if week == 0
}
foreach var of varlist util_tot_pop_* {
	replace `var' = util_tot_pop_nol[2] if week == 0
}

* generate variables under no pandemic
gen byte alive_nop = 1
label var alive_nop "Share alive with no pandemic"
gen byte lockdown_nop = 0
label var lockdown_nop "Lockdown policy with no pandemic"
gen byte s_nop = 1
label var s_nop "Share susceptible (S) with no pandemic"
gen byte i_nop = 0
label var i_nop "Share infected (I) with no pandemic"
gen byte r_nop = 0
label var r_nop "Share recovered (R) with no pandemic"
gen byte d_nop = 0
label var d_nop "Share dead (D) with no pandemic"
sort week
gen float gdp_nop = gdp_nol[1]
label var gdp_nop "GDP with no pandemic"
gen float inv_nop = inv_nol[1]
label var inv_nop "Investment with no pandemic"
gen float util_cons_nop = util_cons_nol[1]
label var util_cons_nop "Flow utility from consumption with no pandemic"
gen float util_life_nop = util_life_nol[1]
label var util_life_nop "Flow utility from being alive with no pandemic"
gen float util_cons_pop_nop = util_cons_pop_nol[1]
label var util_cons_pop_nop "Population flow utility from consumption with no pandemic"
gen float util_life_pop_nop = util_life_pop_nol[1]
label var util_life_pop_nop "Population flow utility from being alive with no pandemic"
gen float util_tot_pop_nop = util_tot_pop_nol[1]
label var util_tot_pop_nop "Total population flow utility with no pandemic"

* normalize GDP and investment to = 1 without pandemic
sort week
foreach var of varlist gdp* {
	sum `var' if week == 0, meanonly
	replace `var' = `var'/gdp_nop
}
foreach var of varlist inv* {
	sum `var' if week == 0, meanonly
	replace `var' = `var'/inv_nop
}

* create "infinite" horizon
set obs 12500
replace week = _n - 1 if week == .

* create cumulative weekly discount factor
global beta_weekly = 1/((1 + 0.05)^(1/52))
gen double beta_weekly_cumulative = ${beta_weekly}^week
label var beta_weekly_cumulative "Cumulative weekly discount factor"

* carry forward values of all variables -- GDP (y_t), investment (x_t), consumption utility (u(c_t)), life utility (nu), share susceptible (S_t), share infected (I_t), share recovered (R_t), share dead (D_t), share alive (1 - S_t - I_t - R_t - D_t), etc.
sum week if lockdown_nol < ., meanonly
local week_max = r(max)
foreach var of varlist * {
	sum `var' if week == `week_max', meanonly
	replace `var' = r(mean) if `var' == .
}

* compute product of total flow utility and cumulative weekly discount factor
egen double util_tot_pop_nop_disc = total(util_tot_pop_nop*beta_weekly_cumulative)
label var util_tot_pop_nop_disc "Discounted sum of utility with no pandemic"
egen double util_tot_pop_nol_disc = total(util_tot_pop_nol*beta_weekly_cumulative)
label var util_tot_pop_nol_disc "Discounted sum of utility with no lockdown"
egen double util_tot_pop_commit_disc = total(util_tot_pop_commit*beta_weekly_cumulative)
label var util_tot_pop_commit_disc "Discounted sum of utility with commitment"
egen double util_tot_pop_nocommit_disc = total(util_tot_pop_nocommit*beta_weekly_cumulative)
label var util_tot_pop_nocommit_disc "Discounted sum of utility with no commitment"

* convert SIRD states to percent
foreach sird in s i r d {
	foreach var of varlist `sird'_nop `sird'_nol `sird'_commit `sird'_nocommit {
		replace `var' = `var'*100
	}
}

* save
save "${DIR_SIMULATIONS}/sim_time_series.dta", replace


*** plot figures
* plot time series of lockdown policy
tw ///
	(connected lockdown_nop lockdown_nol lockdown_commit lockdown_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.2).8, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Lockdown share, L{subscript:t}") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_lockdown, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_lockdown.${fig_format}", replace

* plot time series of GDP
tw ///
	(connected gdp_nop gdp_nol gdp_commit gdp_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(.2(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("GDP, y{subscript:t} (rel. to no pandemic)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_gdp, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_gdp.${fig_format}", replace

* plot time series of aggregate consumption
tw ///
	(connected gdp_nop gdp_nol gdp_commit gdp_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(.2(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Aggregate consumption, c{subscript:t} (rel. to no pandemic)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_cons, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_cons.${fig_format}", replace

* plot time series of investment
tw ///
	(connected inv_nop inv_nol inv_commit inv_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(.2(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Investment, I{subscript:t} (rel. to no pandemic)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_inv, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_inv.${fig_format}", replace

* plot time series of share alive
tw ///
	(connected alive_nop alive_nol alive_commit alive_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(.95(.01)1, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Share alive, 1 - D{subscript:t}") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(8)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_alive, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_alive.${fig_format}", replace

* plot time series of share susceptible, S
tw ///
	(connected s_nop s_nol s_commit s_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(40(10)100, grid gstyle(dot) gmin gmax format(%3.0f)) ///
	xtitle("Week, t") ytitle("Share susceptible, S{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_s, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_s.${fig_format}", replace

* plot time series of share infected, I
tw ///
	(connected i_nop i_nol i_commit i_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(2)10, grid gstyle(dot) gmin gmax format(%2.0f)) ///
	xtitle("Week, t") ytitle("Share infected, I{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(1)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_i, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_i.${fig_format}", replace

* plot time series of share recovered, R
tw ///
	(connected r_nop r_nol r_commit r_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(10)60, grid gstyle(dot) gmin gmax format(%2.0f)) ///
	xtitle("Week, t") ytitle("Share recovered, R{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_r, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_r.${fig_format}", replace

* plot time series of share dead, D
tw ///
	(connected d_nop d_nol d_commit d_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(1)4, grid gstyle(dot) gmin gmax format(%1.0f)) ///
	xtitle("Week, t") ytitle("Share dead, D{subscript:t} (%)") ///
	legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_d, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_d.${fig_format}", replace

* plot time series of flow utility from consumption
tw ///
	(connected util_cons_pop_nop util_cons_pop_nol util_cons_pop_commit util_cons_pop_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Flow utility from consumption, u(c{subscript:t})") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(1)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_u_cons, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_u_cons.${fig_format}", replace

* plot time series of flow utility from being alive
tw ///
	(connected util_life_pop_nop util_life_pop_nol util_life_pop_commit util_life_pop_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Flow utility from being alive, {&nu}(c{subscript:t})") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(1)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_u_life, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_u_life.${fig_format}", replace

* plot time series of total population flow utility
tw ///
	(connected util_tot_pop_nop util_tot_pop_nol util_tot_pop_commit util_tot_pop_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(black blue red green) mcolor(black blue red green) msymbol(i O D T) lpattern(l longdash dash shortdash_dot) lwidth(vthick thick thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Total flow utility, u(c{subscript:t}) + {&nu}") ///
	legend(off) /// legend(order(1 "No pandemic" 2 "No lockdown" 3 "Lockdown with commitment" 4 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_u_tot, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_u_tot.${fig_format}", replace

* plot time series of lockdown policy, commitment vs. no commitment
tw ///
	(connected lockdown_commit lockdown_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Lockdown share, L{subscript:t}") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_lockdown_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_lockdown_c_vs_n.${fig_format}", replace

* plot time series of GDP, commitment vs. no commitment
tw ///
	(connected gdp_commit gdp_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(.2(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("GDP, y{subscript:t} (rel. to no pandemic)") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_gdp_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_gdp_c_vs_n.${fig_format}", replace

* plot time series of aggregate consumption, commitment vs. no commitment
tw ///
	(connected gdp_commit gdp_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.2)1, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Agg. consumption, c{subscript:t} (rel. to no pandemic)") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(5)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_cons_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_cons_c_vs_n.${fig_format}", replace

* plot time series of share susceptible, S, commitment vs. no commitment
tw ///
	(connected s_commit s_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(98(.5)100, grid gstyle(dot) gmin gmax format(%4.1f)) ///
	xtitle("Week, t") ytitle("Share of susceptible, S{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(7)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_s_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_s_c_vs_n.${fig_format}", replace

* plot time series of share infected, I, commitment vs. no commitment
tw ///
	(connected i_commit i_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.05).3, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Share of infected, I{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(1)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_i_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_i_c_vs_n.${fig_format}", replace

* plot time series of share recovered, R, commitment vs. no commitment
tw ///
	(connected r_commit r_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.2)1.2, grid gstyle(dot) gmin gmax format(%2.1f)) ///
	xtitle("Week, t") ytitle("Share of recovered, R{subscript:t} (%)") ///
	legend(off) /// legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_r_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_r_c_vs_n.${fig_format}", replace

* plot time series of share dead, D, commitment vs. no commitment
tw ///
	(connected d_commit d_nocommit week if inrange(week, 0, ${week_plot_max}), lcolor(red green) mcolor(red green) msymbol(D T) lpattern(dash shortdash_dot) lwidth(thick thick) xline(${period_miracle_cure}, lcolor(black) lpattern(.) lwidth(vthick))) ///
	, xlabel(0(10)${week_plot_max}, grid gstyle(dot) gmin gmax format(%2.0f)) ylabel(0(.02).12, grid gstyle(dot) gmin gmax format(%3.2f)) ///
	xtitle("Week, t") ytitle("Share of dead, D{subscript:t} (%)") ///
	legend(order(1 "Lockdown with commitment" 2 "Lockdown without commitment") symxsize(*.4) region(fcolor(none) lcolor(none)) cols(1) ring(0) position(11)) ///
	scale(1.33) plotregion(lcolor(black) margin(l=0 r=0 b=0 t=0)) graphregion(color(white)) /// ysize(2.8) xsize(7) scale(1.9)
	name(sim_d_c_vs_n, replace)
if $fig_export graph export "${DIR_FIGURES}/sim_d_c_vs_n.${fig_format}", replace
