/*--------------------------------------------------------------------
  01_mmrm_primary.sas -- Aim 1 double-programming of R/03_mmrm_primary.R
  MMRM: CHG = TX*AVISIT + BASE + AGEYR + EDCCNTU + APOE4_AN
        unstructured covariance, Kenward-Roger df, REML
  Output to reconcile with output/tables/aim1_mmrm_tx_difference.csv
--------------------------------------------------------------------*/
%let root = /home/&sysuserid./a4-longitudinal-mmrm;
libname a4 "&root./data/derived";

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
run;

proc print data=tx_diff noobs label; title "SAS PROC MIXED: Solanezumab - Placebo, change from baseline PACC"; run;

proc export data=tx_diff outfile="&root./output/tables/aim1_sas_tx_difference.csv" dbms=csv replace; run;
proc export data=lsm     outfile="&root./output/tables/aim1_sas_lsmeans.csv"       dbms=csv replace; run;
