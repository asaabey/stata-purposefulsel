//ModelBuilder v1.3

capture macro drop mb_*
global modelindex = 1
mat modelindexmat = (1,1,1,1,1)

global modelindex = 1
global modeldomainindex =10


capture program drop modelbuild
program modelbuild,rclass
	version 15.1
	args varlist depvar method
	tempname retuvmat retmvmat retuimat retmimat
	
	global primaryvarlist "`varlist'"
	//Define sample
	regmb "`varlist'" "`depvar'" 0.05 `method' multi show defsample
	global modeldomainindex = $modeldomainindex + 1
	
	//Univariate
	di in blue "-----------------------------------------"
	di in blue "|   Uni variate regression               |"
	di in blue "-----------------------------------------"
	local uvvarlist `varlist'
	
	regmb "`varlist'" "`depvar'" 0.2 `method' uni show sample1
	global modeldomainindex = $modeldomainindex + 1
	
	mat `retuvmat'=r(retmat)
	local uvlistout=r(pipevarlist)
	
	//Multivariate
	di in blue "-----------------------------------------"
	di in blue "|   Multi variate regression               |"
	di in blue "-----------------------------------------"
	regmb "`varlist'" "`depvar'" 0.05 `method' multi show sample1
	global modeldomainindex = $modeldomainindex + 1
	mat `retmvmat'=r(retmat)
	local mvlistout=r(pipevarlist)
	
	//Build preliminary main effects
	di in blue "-----------------------------------------"
	di in blue "|   preliminary main effects             |"
	di in blue "-----------------------------------------"
	BuildPrelimModel,rmat(r(retmat))
	local varlist5 `r(pipevarlist)' 
	global modeldomainindex = $modeldomainindex + 1
	local mvprelimlistout=r(pipevarlist)
	
	
	//Multivariate-uni interaction
	di in blue "-----------------------------------------"
	di in blue "|   interaction only                     |"
	di in blue "-----------------------------------------"
	regmb_int_only "`r(pipevarlist)'" "`depvar'" 0.05 `method' multi_int hide
	local varlist6 "`varlist5' `r(pipevarlist)'"
	di "pre-final model :`varlist6'"
	global modeldomainindex = $modeldomainindex + 1
	mat `retuimat'=r(retmat)
	local mvvarlist=r(pipevarlist)
	local mivarlist "`mvprelimlistout' `mvvarlist'"
	
	
	//Multivariate-multi interaction
	di in blue "-----------------------------------------"
	di in blue "|   interaction & main effects           |"
	di in blue "-----------------------------------------"
	regmb "`varlist6'" "`depvar'" 0.05 `method' multi show sample1
	di "Final model : `r(pipevarlist)'"
	mat `retmimat'=r(retmat)
	local fvarlist=r(pipevarlist)
	
	//rclass returns
	return mat univarmat=`retuvmat'
	return mat multivarmat=`retmvmat'
	return mat uniintmat=`retuimat'
	return mat multiintmat=`retmimat'
	
	return local univarlist "`uvvarlist'"
	return local multivarlist "`mvvarlist'"
	return local uniintlist "`uivarlist'"
	return local multiintlist "`fvarlist'"
end


capture program drop regmb_int_only
program regmb_int_only,rclass
	version 15.1
	args varlist depvar p_thresh method type hide_non_sig 
	tempname retmat 
	mat `retmat'=(1,1,1,1,1)

	
	foreach va of local varlist {
		foreach vb of local varlist {
			if `va'!=`vb' {
				
				resolveWeights,primaryvarlist($primaryvarlist) var("`vb'")
				local vbi=r(retvar)
				resolveWeights,primaryvarlist($primaryvarlist) var("`va'")
				local vai=r(retvar)
				//local int_term="`vb'#`va'"
				local int_term="`vbi'#`vai'"
				regmb "`varlist' `int_term'" `depvar' `p_thresh' `method' `type' `hide_non_sig' sample1
				local pvl=r(pipevarlist)
				mat rt1=r(retmat)
				local rowsmat=rowsof(rt1)
				mat `retmat'=(`retmat'\rt1[`rowsmat'..`rowsmat',1..5])
				
			}
			if "`pvl'"!="." {
					local varlist3 "`varlist3' `pvl'"
			}
			
		}
	}
	
	mat `retmat'=`retmat'[2...,1...]
	mat colnames `retmat'=coef llcoef ulcoef p sigflag
	return matrix retmat=`retmat'
	return local pipevarlist `varlist3'
end

capture program drop regmb
program regmb,rclass 
	version 15.1
	args varlist depvar p_thresh method type hide_non_sig defsample
	local disptemplate="`type' : Covariate : `var', p : `p' `sigflag'"
	tempname retmat 
	mat `retmat'=(1,1,1,1,1,1)
	local pipevarlist ""
	local rowvarlist ""
	local defsamplestr " if e(sample)"
	local sigflag ""
	local nsigflag 0
	local keepflag 0
	local i=1
	local enumv=1
	local ell
	
	if "`defsample'"=="defsample" {
		local defsamplestr ""
	}
	if "`method'"=="stcox" {
			local depvar ""
	}
	if "`type'"=="multi" | "`type'"=="multi_int" {
	

		
		
		local varlist=subinstr("`varlist'","[k]","",.)
		
		local varlist=subinstr("`varlist'"," .","",.)
		
		quietly  `method' `depvar' `varlist' `defsamplestr'
		//di in red "============ ll ->" e(ll) " ====== df :" e(df_m) " ==== N "  e(N)
		tempname rtbl 
		
		mat `rtbl'=r(table)
		mat `rtbl'=`rtbl''
		local ell:di e(ll)
		
		
	}
	foreach var of local varlist {
		
		
		//check for keep flag [k]
		if strpos("`var'","[k]")>0 {
			local keepflag 1
			//strip the keepflag
			local var=subinstr("`var'","[k]","",1)	
		}
		
		if "`type'"=="uni" | missing("`type'"){
			
			quietly `method' `depvar' `var' `defsamplestr'
			tempname rtbl 
		
			mat `rtbl'=r(table)
			mat `rtbl'=`rtbl''
			local i=1
			local ell:di e(ll)
			
		}
		else if "`type'"=="multi_int" {
			if "`method'"=="stcox" {
					local i=rowsof(`rtbl')-0
			}
			else {
					local i=rowsof(`rtbl')-1
			}
			
		}
		
		local p:di %5.2f `rtbl'[`i',4]
		local coef:di %5.2f `rtbl'[`i',1]
		local llcoef:di %5.2f `rtbl'[`i',5]
		local ulcoef:di %5.2f `rtbl'[`i',6]
		
		local rowvarlist `rowvarlist' `var'
		
		if `p'< `p_thresh' & `keepflag'==0{
			local sigflag "***"	
			local nsigflag 1
			//di in red "DIAG::::: `type'"
			if "`type'"=="multi_int" {
				
				//di in red "DIAG::::: `var'"
				if strpos("`var'","#")>0 {
					local pipevarlist `pipevarlist' `var'
				}
			} 
			else {
				local pipevarlist `pipevarlist' `var'
			}
			
		} 
		else {
			local sigflag ""
			local nsigflag 0
		}
		
		mat b= (`rtbl'[`i',1],`rtbl'[`i',5],`rtbl'[`i',6],`rtbl'[`i',4],`nsigflag',`enumv')
		mat `retmat'= (`retmat'\b)
		
		// Console output : Start ----------------------
		if "`type'"!="multi_int" {
			if "`sigflag'"=="***" {
					stdout `method' `type' `var' `coef' `llcoef' `ulcoef' `p' `sigflag'
			}
			else {
					if "`hide_non_sig'"!="hide" {
						stdout `method' `type' `var' `coef' `llcoef' `ulcoef' `p' `sigflag'
					}
			}
			
		} 
		else if "`type'"=="multi_int"{
			if strpos("`var'","#")>0 {
				if "`sigflag'"=="***" {
					stdout `method' `type' `var' `coef' `llcoef' `ulcoef' `p' `sigflag'

				}
				else {
					if "`hide_non_sig'"!="hide" {
						stdout `method' `type' `var' `coef' `llcoef' `ulcoef' `p' `sigflag'
					}
				}
			}
		}
		// Console output : End ----------------------
		local ++i
		local ++enumv
		mat drop b
		local drop p
		
	}
	mat `retmat'=`retmat'[2...,1...]
	mat colnames `retmat'=coef llcoef ulcoef p sigflag enumv
	mat rownames `retmat'=`rowvarlist'
	
	if "`type'"!="multi_int" {
		registerModel,modelvarlist("`pipevarlist'") method("`method'") prefix(12) df(`e(df_m)') ll(`e(ll)')  ll0(`e(ll_0)') n(`e(N)')
	}
	
	 
	
	return local df=e(df_m)
	return local ll=e(ll)
	return local ll_0=e(ll_0)
	
	return matrix retmat=`retmat'
	return local pipevarlist `pipevarlist'
	
	
	
	
end

capture program drop stdout
program stdout
	version 15.1
	args method type covariate coef llcoef ulcoef p sigflag
	local txthl ""
	if "`sigflag'"=="***"{
		local txthl "in yellow"
	}
	//display `txthl' "method:`method',type:`type' ,covariate:`covariate' ,coef:`coef' [`llcoef'-`ulcoef'] p:`p' ,`sigflag'"
end

capture program drop GetVarName
program GetVarName,rclass
	version 15.1
	syntax,varlist(string) enum(int)
	local i=1
	foreach v of local varlist {
		if `i'==`enum' {
			return local retvar `v'
		}
		local ++i
	}
end


capture program drop BuildPrelimModel
program BuildPrelimModel, rclass
	version 15.1
	syntax,rmat(string)
	
	local depvar "rrt"
	local method "logit"
	local deltab_thresh=0.2
	tempname mf msb msfp mssp
	mat `mf'=`rmat'
	local mfrows = rowsof(`mf')
	local varlista: rownames(`mf')
	mat `msb'=(1,1,1,1,1,1)
	local vla ""
	local va ""
	forval i=1/`mfrows' {
		if `mf'[`i',5]==0 {
			mat `msb'=`msb' \ `mf'[`i'..`i',1...]
		}
	}
	mat `msb'=`msb'[2...,1...]
	mata : st_matrix("`msb'",sort(st_matrix("`msb'"),-4))
	
	mat list `msb'
	
	local msbrows=rowsof(`msb')
	local varlistb ""
	// Generate non-sig list
	forval i=1/`msbrows' {
		local v_ind=`msb'[`i',6]
		GetVarName,varlist("`varlista'") enum(`v_ind')
		local varlistb `varlistb' `r(retvar)'
	}
	
	local varlistc `varlista'
	local fp_varlist `varlista'
	local sp_varlist ""
	
	forval i=1/`msbrows' {
		local v_ind=`msb'[`i',6]
		GetVarName,varlist("`varlista'") enum(`v_ind')
		local varexcl `r(retvar)'
		local vardrop=0
		di "==== EXCLUSION VAR : `varexcl' ======"
		
		
		local varlistc:di strtrim(subinstr("`varlistc'","`varexcl'","",1))
		local fp_varlist_temp `fp_varlist'
		
		
		local sp_varlist:di strtrim(subinstr("`fp_varlist_temp'","`varexcl'","",1))
		//first pass
		di "==== FIRSTPASS `varexcl' ======"
		di "==== FIRSTPASS VARLIST : `fp_varlist' ======"
		regmb "`fp_varlist'" "`depvar'" 0.05 `method' multi show sample1
		mat `msfp'=r(retmat)
		//mat list `msfp'
		
		
		//second pass
		di "==== SECONDPASS `varexcl' ======"
		di "==== SECONDPASS VARLIST : `sp_varlist' ======"
		regmb "`sp_varlist'" "`depvar'" 0.05 `method' multi show sample1
		mat `mssp'=r(retmat)
		//mat list `mssp'
		
		local msfprows=rowsof(`msfp')
		local mssprows=rowsof(`mssp')
		
		di "==== STARTING LOOP ======"
		forval i=1/`msfprows' {
			local fp_ind=`msfp'[`i',6]
			local fp_sig=`msfp'[`i',5]
			local fp_b=`msfp'[`i',1]
			
			forval j=1/`mssprows' {
					local sp_ind=`mssp'[`j',6]
					local sp_sig=`mssp'[`j',5]
					local sp_b=`mssp'[`j',1]
					
					//di "fp_ind: `fp_ind'; sp_ind: `sp_ind'; fp_sig: `fp_sig'"
					
					GetVarName,varlist("`fp_varlist'") enum(`fp_ind')
					local varname_fp `r(retvar)' 
					
					GetVarName,varlist("`sp_varlist'") enum(`sp_ind')
					local varname_sp `r(retvar)' 
					
					//di "--------- `varname_fp' vs. `varname_sp' : `fp_sig'----------------"
					// ERROR
					//if `fp_ind'==`sp_ind' & `fp_sig'==1 {
					//if "varname_fp"=="varname_sp" & `fp_sig'==1 {
					if "`varname_fp'"=="`varname_sp'" & `fp_sig'==1 {
						local deltab = abs((`fp_b'-`sp_b')/(`fp_b'))
						
						//GetVarName,varlist("`varlista'") enum(`sp_ind')
						//GetVarName,varlist("`sp_varlist'") enum(`sp_ind')
						//local varfocus `r(retvar)'
						di "--------- `varname_fp'(`fp_ind') vs. `varname_sp'(`sp_ind') ----------------"
						di "var fp : `i'; var sp : `j'"
						di "fp_b:`fp_b'; sp_b:`sp_b'"
						di "----db:  `deltab'"
						di "fp: `fp_varlist'"
						di "sp: `sp_varlist'"
						
						if `deltab'>`deltab_thresh' {
							// cov cannot be dropped
							//local fp_varlist `fp_varlist'
							local vardrop=0
						}
						else {
							// cov can be dropped
							//local fp_varlist `sp_varlist'
							
							local vardrop=1
						}
						
					}
			}
			
		}
		
		
		if `vardrop'==1 {
			local fp_varlist `sp_varlist'
			di "===`varexcl' dropped === fp: `fp_varlist'"
		}
		else {
			local fp_varlist `fp_varlist'
			di "===`varexcl' NOT dropped ==='fp: `fp_varlist'"
		}
		
		
	}
	
	
	return local vla `varlista'
	return local vlb `varlistb'
	return local vlc `varlistc'
	return local pipevarlist `fp_varlist'
	return matrix msb=`msb'
end

// Register model in the Model Index Matrix (modelindexmat)
// Dependent vars: Globals
// Dependent prog: None
// Called by : regmb
capture program drop registerModel
program registerModel
	version 15.1
	syntax,modelvarlist(string) method(string) prefix(int) df(int) ll(string) ll0(string) n(int)
	tempname tmprow
	
	local i=($modeldomainindex*1000)+$modelindex
	mat `tmprow'=(`i',`df',`ll',`ll0',`n')
	mat modelindexmat = (modelindexmat \ `tmprow')
	
	global mb_lmodel`i' "`modelvarlist'"
	
	global modelindex = $modelindex + 1
end

capture program drop resolveWeights
program resolveWeights,rclass
	syntax,primaryvarlist(string) var(string)
	foreach v of local primaryvarlist {
		if strpos("`v'","`var'")==3 & (strlen("`v'")-strlen("`var'"))==2 {
				return local retvar "`v'"
		}
		else if strpos("`v'","`var'")==1 & (strlen("`v'")-strlen("`var'"))==0 {
				return local retvar "`v'"
		}
	}
	
end

