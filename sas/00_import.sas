/*--------------------------------------------------------------------
  00_import.sas -- import derived analysis datasets exported from R
  Run R/06_export_for_sas.R first; it writes data/derived/*.csv
  Tested on SAS OnDemand for Academics (SAS 9.4M8).
--------------------------------------------------------------------*/
%let root = /home/&sysuserid./a4-longitudinal-mmrm;   /* adjust to your SAS ODA path */
libname a4 "&root./data/derived";

proc import datafile="&root./data/derived/adqs_pacc.csv" out=a4.adqs dbms=csv replace; guessingrows=max; run;
proc import datafile="&root./data/derived/adsl.csv"      out=a4.adsl dbms=csv replace; guessingrows=max; run;

proc contents data=a4.adqs varnum short; run;
