local tblprf "logit"
local tbltxt "Logistic"
local tblnum "3"

//modelbuild "c.egfr c.acrln c.age_c male dm htn c.hb_c c.alb_c c.potassium_c c.bicarb_c zone slopeh sbp" rrt `tblprf'

//modelbuild "c.egfr c.acrln c.age_c male dm htn c.hb_c c.alb_c c.potassium_c c.bicarb_c slopeh sbp" rrt `tblprf'
modelbuild "c.egfr5 c.acrln c.age_c male dm htn c.hb_c c.alb_c c.potassium_c c.bicarb_c slopeh sbp" rrt `tblprf'
tempname m1 m2 m3 m4
mat `m1'=(r(univarmat))
mat `m2'=(r(multivarmat))
mat `m3'=(r(uniintmat))
mat `m4'=(r(multiintmat))


applylabmat, rmat(`m1')
tempname m
mat `m'=(r(ret))
//applymatfilter, rmat(`m') filtercol(5) filtercriteria(1)
//mat `m'=(r(retmat))
putdocx paragraph
putdocx text ("Table `tblnum'.1 `tbltxt' regression, Univariate analysis"),font($fnameNm,$fsizeNm)
putdocx table `tblpfr'univartbl = matrix(`m'),nformat(%9.3f) memtable rownames colnames border(start, nil) border(insideH, nil) border(insideV, nil) border(end, nil) 
moddoctbl , rmat(`m') tblname(`tblpfr'univartbl)
mat drop `m'

applylabmat, rmat(`m2')
tempname m
mat `m'=(r(ret))
applymatfilter, rmat(`m') filtercol(5) filtercriteria(1)
mat `m'=(r(retmat))
putdocx paragraph
putdocx text ("Table `tblnum'.2 `tbltxt' regression, Multivariate analysis")
putdocx table `tblpfr'multivartbl = matrix(`m'),nformat(%9.3f) memtable  rownames colnames note("Note: Table 2.2")
moddoctbl , rmat(`m') tblname(`tblpfr'multivartbl)
mat drop `m'

mat list `m3'
applylabmat, rmat(`m3')
tempname m
mat `m'=(r(ret))
applymatfilter, rmat(`m') filtercol(5) filtercriteria(0)
mat `m'=(r(retmat))
putdocx paragraph
putdocx text ("Table `tblnum'.3 `tbltxt' regression, Interaction analysis")
putdocx table `tblpfr'uniinttbl = matrix(`m'),nformat(%9.3f) memtable rownames colnames border(start, nil) border(insideH, nil) border(insideV, nil) border(end, nil)  //note("Note: Table 2.1")
moddoctbl , rmat(`m') tblname(`tblpfr'uniinttbl)
mat drop `m'


applylabmat, rmat(`m4')
tempname m
mat `m'=(r(ret))
applymatfilter, rmat(`m') filtercol(5) filtercriteria(1)
mat `m'=(r(retmat))
putdocx paragraph
putdocx text ("Table `tblnum'.4 `tbltxt' regression, Final Multivariate model with interactions")
putdocx table `tblpfr'multiinttbl = matrix(`m'),nformat(%9.3f) memtable rownames colnames note("Note: Table 2.3")
//moddoctbl `tblpfr'multiinttbl
moddoctbl , rmat(`m') tblname(`tblpfr'multiinttbl)
mat drop `m'

logit
est store final_lm


