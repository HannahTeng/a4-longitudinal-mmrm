/*--------------------------------------------------------------------
  01_mmrm_primary.sas -- Aim 1 double-programming of R/03_mmrm_primary.R
  MMRM: CHG = TX*AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN
        unstructured covariance, Kenward-Roger df, REML
  Output to reconcile with output/tables/aim1_mmrm_tx_difference.csv
--------------------------------------------------------------------*/
%let root = /home/&sysuserid./A4;   /* flat SAS ODA folder holding the CSVs and these scripts */
libname a4 "&root.";

data mitt;
  set a4.adqs;
  where SUBSTUDY = "A4" and MITTFL = 1 and AVISITN > 0
        and not missing(CHG) and not missing(BASE)
        and not missing(AGEYR) and not missing(EDCCNTU) and not missing(APOE4_AN);
run;

ods output LSMeans=lsm Diffs=diffs CovParms=covp SolutionF=solf;
proc mixed data=mitt method=reml;
  class BID TX(ref="Placebo") AVISIT APOE4_AN(ref="non-carrier");
  model CHG = TX AVISIT TX*AVISIT BASE AGEYR EDCCNTU APOE4_AN / solution ddfm=kr;
  repeated AVISIT / subject=BID type=un;
  lsmeans TX*AVISIT / diff cl alpha=0.05;
run;

/* keep only within-visit treatment contrasts (Solanezumab - Placebo) */
data tx_diff;
  set diffs;
  where AVISIT = _AVISIT and TX = "Solanezumab" and _TX = "Placebo";
  keep AVISIT Estimate StdErr DF Lower Upper Probt;
  format Estimate StdErr DF Lower Upper best16. Probt best16.;   /* full precision for CSV export (default d8.4 rounds to 4 dp) */
run;

proc print data=tx_diff noobs label; title "SAS PROC MIXED: Solanezumab - Placebo, change from baseline PACC"; run;

/* full-precision CSV export; data _null_ + put ignores the d8.4 formats that PROC EXPORT would inherit */
data _null_;
  set tx_diff;
  file "&root./aim1_sas_tx_difference.csv" dsd lrecl=32767;
  if _n_ = 1 then put "AVISIT,Estimate,StdErr,DF,Probt,Lower,Upper";
  put AVISIT :$8. Estimate :best32. StdErr :best32. DF :best32. Probt :best32. Lower :best32. Upper :best32.;
run;

data _null_;
  set lsm;
  file "&root./aim1_sas_lsmeans.csv" dsd lrecl=32767;
  if _n_ = 1 then put "Effect,TX,AVISIT,Estimate,StdErr,DF,tValue,Probt,Alpha,Lower,Upper";
  put Effect :$20. TX :$20. AVISIT :$8. Estimate :best32. StdErr :best32. DF :best32. tValue :best32. Probt :best32. Alpha :best32. Lower :best32. Upper :best32.;
run;
