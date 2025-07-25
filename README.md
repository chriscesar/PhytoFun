# PhytoFun
Examinations of phytoplankton time series data

Incorporating phytoplankton abundances and carbon contents.

## TO DO

* ~Check conversions in zooplankton and phytoplankton abundances~
* Check 'missing' taxa for carbon values.  Can we use congeners/parent taxa values?
* Clean code into a single script

## Workflow
* script: '00Newsetup.R'
    + Read Phyto data since 2000 from EA BOXI extract & save as a smaller .Rdata object

* script: '01_DataImportandFormat_v2.R'
    + Extract unique taxon names from the phyto data and use worrms::wm_name2id to assign Aphia ID values to taxa.
  Remove NA values (unassigned taxa) and map WORMS record info to each Aphia match.
  Export 'NA' Aphia values for manual match ups.
  Combine both records into a single 'records' object.
  Left join the taxon record info to the phytoplankton abundance data.
  Write out objects.
    + Import carbon-by-taxon data.  Based on PML & EA estimates. Where a taxon has >1 value, mean and median values are calculated.
  Carbon values appended to phyto data and abundance values multiplied by carbon content per cell values.
    + Summary of taxa with 'missing' carbon values exported.
    + Data exported.
