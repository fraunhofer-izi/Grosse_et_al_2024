# Repository Layout

```         
assets/
├── manifest.yaml                       <------ Analysis-wide settings and paths (template, insert your own paths!)
├── manifest_full.yaml                  <------ Same as the above, with samples that were not included in the publication present
├── README.md                           <------ This file
├── cellCycleMarkers.R                  <------ Cycle markers for basic QC and ...
├── houseKeepingMarkers.R               <------ Housekeeping markers for basic QC from Tirosh et al. (2016), doi:10.1126/science.aad0501
├── annotations/                        <------ Histological annotations transferred via Loupe browser
├── IF_pictures_for_merged_panels/      <------ IF picture for incorporation into reports
├── single_cell_referece/               <------ Script to download and cell type annotation for Stubenvoll et al. (2025) data
├── GTEx_v9_snRNA_seq_Eraslan_et_al_2022/               <------ GTEx v9 data from Eraslan et al. (2022), subset to skin and converted for convenience
└── TCIA-CellTypeFractionsData_SKCM_only_march2024.tsv  <------ TCIA cell type fraction data for TCGA cohort from https://tcia.at/cellTypeFractions (Charoentong et al. 2016, doi:10.1016/j.celrep.2016.12.019)
```
