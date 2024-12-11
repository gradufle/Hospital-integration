clear
global in C:\Users\lguerin\OneDrive - ehesp.fr\Dossier partage\PaperGHT\FINALWPGHT //Setup path for datasource
use "${in}/datasetup3.dta", clear

****Setup*****

global dep c_iftr_chain 												//Dependant variable
global tit y_hat te p_r80 p_charlson_glob p_age c_ci q_rh30s p_bedblock // Covariates
global date 2016 														//Policy implementation date-1
global tip optimized projected 											// type of estimation method (Arkhangelsky / Kranz)
global bot 2  															//number of bootstrap repetitions
global cut 2020 														// cutoff for sample: max=2023; min=2014 can imput multiple cutoffs

tab oldght, m
drop if oldght==4 | oldght==5 											//Keep only relevent hospitals. If you want to drop private add | oldght==6

xtset id i_year

******TWFE Estimates of the Determinants of the IFTR, 2013-2016 (Table 2):******

foreach v in $dep{	
	xtreg `v' $tit i.i_year if i_year<=2016, fe robust
	est sto R
foreach i in 1 2 3 6 7{
			xtreg `v' $tit i.i_year if oldght==`i' & i_year<=2016, fe robust
			est sto R`i'
}
}

	# delimit ;	
	esttab R R* using "${in}/regcov2.csv", replace plain cells(b(star fmt(%9.3f))se(fmt(%9.3f))) noomitted nobase lab compress  
	varwidth(40) stats(N r2_o r2_w, fmt(%9.0f %9.3f %9.3f) labels("Observations" "R2 overall" "R2 within")) ;# delimit cr

******Synthetic Difference in differences estimation loop (Table 3,4,5):******

*****Foreach cutoff point:*****
foreach p in $cut{
****Foreach estimation method:****
foreach n in $tip{
***Foreach dependant variable:***
foreach v in $dep{
	
**For Total sample:**
    preserve
	//Fully balanced panel:	
	qui xtreg `v' $tit, fe 
	capture gene sample2 =e(sample)
	drop if sample2 != 1
	keep if i_year<`p'
	//Create treated/control group:
	capture drop streated
	gen streated = 0
	replace streated=1 if pub==1 & i_year>$date & oldght!=7
	//Fully balanced panel:
	bysort id (i_year): gene nb = _N
		qui sum nb
		keep if nb==`r(max)'
		drop nb 		
	disp "`i'" "`v'" "`n'" "`p'"
*Model estimation:*
sdid `v' id i_year streated, method(sdid) vce(bootstrap) reps($bot) covariates($tit, `n') graph
	est sto j`i'`n'`p'tot	
	restore
	
**For Heterogeneity analysis:**	
	forvalues i=1(1)3{
	preserve
	//Fully balanced panel:
	qui xtreg `v' $tit, fe 
	capture gene sample2 =e(sample)
	drop if sample2 != 1
	keep if i_year<`p'
	//Create treated/control group:
	capture drop streated
	gen streated = 0
	replace streated=1 if oldght==`i' & i_year>$date
	drop if pub==1 & oldght!=`i' & oldght!=7 		
	sum streated 
	//Fully balanced panel:
	bysort id (i_year): gene nb = _N
		qui sum nb
		keep if nb==`r(max)'
		drop nb 	
	
	disp "`i'" "`v'" "`n'" "`p'"
*Model estimation:*
sdid `v' id i_year streated, method(sdid) vce(bootstrap) reps($bot) covariates($tit, `n') graph g1_opt(ylabel(-20(1)20) scheme(plottig)) g1on g2_opt(ylabel(0(1)4) legend(position(3)cols(1)) scheme(plottig))	
	est sto j`i'`n'`p'	
	restore
	}
	}
	}
	}
	
		# delimit ;	
		esttab j* using "${in}/resj2.csv", replace plain cells(b(star fmt(%9.3f))se(fmt(%9.3f))) noomitted nobase lab compress  
		varwidth(40) stats(N r2_o r2_w, fmt(%9.0f %9.3f %9.3f) labels("Observations" "R2 overall" "R2 within")); # delimit cr	
	