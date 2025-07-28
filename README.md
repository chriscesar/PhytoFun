# PhytoFun
Examinations of phytoplankton time series data

Incorporating phytoplankton abundances and carbon contents.

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

* script: '02_AppendLifeforms.R'
    + Loads phyto data generated in '01_DataImportandFormat_v2.R' and correctly-formats dates
    + Extracts BIOSYS site code from site_station_name variable (currently Site name with BIOSYS code appended)
    + Loads lifeform data from PLET Master (v7), amended to remove missing values
    + Generate and export list of taxa with their Lifeform values & Carbon values

* script: '03_JoinPhytoZoops.R'
    + Loads zooplankton data generated in 'https://github.com/chriscesar/ZoopData_Expl' and trims BIOSYS codes
    + Loads phyto data generated in '02_AppendLifeforms.R' & trims it to retain the BIOSYS sites included in the zooplankton data
    + Trims variables across both phyto & zoop data, homogenises variable names and adds 'units' variables for clarity
    + Binds phyto and zoop data (long format) into a singe df & appends site metadata by left_join of zooplankton metadata
    + Exports


## TO DO

* ~Check conversions in zooplankton and phytoplankton abundances~
* ~Check 'missing' taxa for carbon values.  Can we use congeners/parent taxa values?~
* ~Update script '01_....' to incorporate amended carbon values~
* ~Update/sense check subsequent scripts to ensure no errors~
* Apparent under-estimation of phyto carbon (or OVER-estimation of zoops)
* ?Clean code into a single script?