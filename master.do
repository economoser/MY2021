********************************************************************************
* MASTER  Create time-series and comparative-statics figures for an economy
* subject to optimal lockdown policy with and without commitment for Moser and
* Yared (2021). Stata master file that can also call master.m via MATLAB shell
* -- run this second!
********************************************************************************


*** initial housekeeping
* settings
set more off
macro drop _all
clear all
timer clear 1
timer on 1
set seed 1
set type double
set excelxlsxlargefile on
set graphics on
set varabbrev off
set rmsg on
set matsize 11000
set linesize 100
cap log close _all

* directories
foreach dir in ///
	"SET_MAIN_DIRECTORY_HERE" ///
	{
	cap confirm file "`dir'"
	if !_rc global DIR_MAIN = "`dir'"
}
if "${DIR_MAIN}" == "" {
	disp as error "USER ERROR: Must specify main directory in line 28!"
	error 1
}
global DIR_CODE = "${DIR_MAIN}/_code"
global DIR_LOG = "${DIR_MAIN}/_log"
global DIR_RESULTS = "${DIR_MAIN}/_results"
global DIR_FIGURES = "${DIR_RESULTS}/_figures"
global DIR_SIMULATIONS = "${DIR_RESULTS}/_simulations"
foreach dir in ///
	"${DIR_CODE}" ///
	"${DIR_LOG}" ///
	"${DIR_RESULTS}" ///
	"${DIR_FIGURES}" ///
	"${DIR_SIMULATIONS}" ///
	{
	cap confirm file "`dir'"
	if _rc {
		!mkdir "`dir'"
	}
}
foreach dir in ///
	"/Applications/MATLAB_R2020b.app/bin" ///
	"/Applications/MATLAB_R2021a.app/bin" ///
	{
	cap confirm file "`dir'"
	if !_rc global DIR_MATLAB = "`dir'/"
}

* switches
global run_simulations = 1 // run simulations in MATLAB
global run_figures_time_series = 1 // create figures for time series simulations in Stata
global run_figures_comp_stat = 1 // create figures for comparative statics in Stata
global fig_export = 1 // 0 = do not export figures, 1 = export figures
global fig_format = "eps" // format of figures to be exported ("eps" or "pdf")

* open log file
log using "${DIR_LOG}/log_covid_stata.log", text name(log_covid) replace


*** run MATLAB code
if $run_simulations {
	!${DIR_MATLAB}matlab -nodesktop -nodisplay <"${DIR_CODE}/master.m"
}


*** create time-series figures
if $run_figures_time_series do "${DIR_CODE}/time_series.do"


*** create comparative-statics figures
if $run_figures_comp_stat do "${DIR_CODE}/comp_stat.do"


*** final housekeeping
timer off 1
timer list 1
disp "FINISHED ON ${S_DATE} AT ${S_TIME} IN A TOTAL OF `=r(t1)' SECONDS."
log close _all
