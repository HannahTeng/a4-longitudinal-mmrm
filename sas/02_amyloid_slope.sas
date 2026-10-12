/*--------------------------------------------------------------------
  02_amyloid_slope.sas -- Aim 2 double-programming of R/04_amyloid_slope.R
  (a) random intercept/slope model on PACC composite (A4 placebo + LEARN)
  (b) covariance-structure comparison for the MMRM on change from baseline
--------------------------------------------------------------------*/
%let root = /home/&sysuserid./a4-longitudinal-mmrm;
libname a4 "&root./data/derived";

data obs;
  set a4.adqs;
  where COHORT in ("A4 Placebo", "LEARN (Abeta-)") and not missing(AVAL) and not missing(AMYLCENT)
        and not missing(AGEYR) and not missing(EDCCNTU) and not missing(APOE4_AN);
  YEARS  = AVISITN / 52.18;
  CENT10 = AMYLCENT / 10;
  length AMYL_GRP $8;
  if AMYLCENT < 20 then AMYL_GRP = "<20 CL";
  else if AMYLCENT <= 50 then AMYL_GRP = "20-50 CL";
  else AMYL_GRP = ">50 CL";
run;

/* (a) random slopes -- compare to output/tables/aim2_lmm_random_slope.csv */
ods output SolutionF=aim2_fixed;
proc mixed data=obs method=reml;
  class BID APOE4_AN(ref="non-carrier");
  model AVAL = CENT10 YEARS CENT10*YEARS AGEYR EDCCNTU APOE4_AN / solution ddfm=kr;
  random intercept YEARS / subject=BID type=un;
run;
proc export data=aim2_fixed outfile="&root./output/tables/aim2_sas_lmm_fixed.csv" dbms=csv replace; run;

/* (b) covariance structures on change from baseline (post-baseline visits) */
data obs_chg; set obs; where AVISITN > 0 and not missing(CHG); run;

%macro covfit(type, label, sastype);
  ods output FitStatistics=fs_&type;
  proc mixed data=obs_chg method=reml;
    class BID AVISIT AMYL_GRP APOE4_AN(ref="non-carrier");
    model CHG = AMYL_GRP AVISIT AMYL_GRP*AVISIT BASE AGEYR EDCCNTU APOE4_AN / ddfm=kr;
    repeated AVISIT / subject=BID type=&sastype;
  run;
  data fs_&type; length structure $32; set fs_&type; structure = "&label"; run;
%mend;
%covfit(cs,   Compound symmetry,   cs);
%covfit(ar1,  AR(1),               ar(1));
%covfit(arh1, AR(1) heterogeneous, arh(1));
%covfit(un,   Unstructured,        un);

data covcomp; set fs_cs fs_ar1 fs_arh1 fs_un; run;
proc print data=covcomp noobs; title "Covariance structure comparison (-2 REML LL, AIC, BIC)"; run;
