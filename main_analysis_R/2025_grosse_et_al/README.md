# Description

# Reproduction

For a full reproduction, follow the steps outlined in `Interactive_analysis.R`. 
Note that you will need to download all input files and enter their locations 
into the `assets/manifest.yaml` file.

If you only want to regenerate the final pictured from existing intermediate data,
enter the path information in the manifest run only the library calls at the top 
and the lines to build report pages 15 an 16 at the bottom of `Interactive_analysis.R`.

Note that in any case, you will likely need to install `renv` and restore all depedencies
when you first open an R session in this project. Usually, this should be no harder 
than calling `renv::restore()` from where the `renv.lock` file is. In case of problems 
please consult the `renv` documantation ()[https://rstudio.github.io/renv/articles/renv.html]
or open an issue in this repo.

# Repository Layout

```         
2025_grosse_et_al/
├── analysis/               <---------- R Markdown files for reports
├── code/                   <---------- R scripts to process data
├── assets/                 <---------- Manifest and auxilliary files that are synchronized
│   ├── manifest.yaml           <------ Analysis-wide settings and paths
│   ├── README.md               <------ See for detailed descriptions of contents
│   └── ...
├── data                    <---------- Caches for various downloaded datasets etc., not synchronized
├── public/                 <---------- Report pages are built here
├── renv.lock               <---------- Dependency version control
├── README.md
├── .gitignore
├── .Rprofile               <---------- Code run on R session startup, initialized version control
├── 2025_grosse_et_al.Rproj <---------- RStudio project marker
└── Interactive_analysis.R  <---------- Read this if you want to redo the analysis
```
