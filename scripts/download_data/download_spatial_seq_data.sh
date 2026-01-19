#!/bin/bash

# download dataset GSE314509 from GEO

mkdir -p ../../data/spatial_seq_data/
cd ../../data/spatial_seq_data/
wget -nc -O "GSE314509_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE314509&format=file"
tar -xvf GSE314509_RAW.tar

# Separate files per sample into separate folders
mkdir HD_PM1 HD_PM2 HD_PM3 HD_PM4 HD_MM1 HD_MM2 HD_MM3
mv GSM*HD_PM1* HD_PM1
mv GSM*HD_PM2* HD_PM2
mv GSM*HD_PM3* HD_PM3
mv GSM*HD_PM4* HD_PM4
mv GSM*HD_MM1* HD_MM1
mv GSM*HD_MM2* HD_MM2
mv GSM*HD_MM3* HD_MM3
