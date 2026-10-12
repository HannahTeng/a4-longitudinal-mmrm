/*--------------------------------------------------------------------
  00_import.sas -- import derived analysis datasets exported from R
  Run R/06_export_for_sas.R first; it writes data/derived/*.csv
  Tested on SAS OnDemand for Academics (SAS 9.4M8).
--------------------------------------------------------------------*/
%let root = /home/&sysuserid./A4;   /* flat SAS ODA folder holding the CSVs and these scripts */
libname a4 "&root.";

proc import datafile="&root./adqs_pacc.csv" out=a4.adqs dbms=csv replace; guessingrows=max; run;
proc import datafile="&root./adsl.csv"      out=a4.adsl dbms=csv replace; guessingrows=max; run;

proc contents data=a4.adqs varnum short; run;
