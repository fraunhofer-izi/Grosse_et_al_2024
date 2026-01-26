#!/bin/bash

# download dataset GSE277165 from geo (supplementary .tar of raw data)

wget -nc -O "GSE277165_RAW.tar" "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE277165&format=file"
mkdir -p ./GSE277165_extracted/
cd GSE277165_extracted
tar -xvf ../GSE277165_RAW.tar
cd ..
wget "https://ftp.ncbi.nlm.nih.gov/geo/series/GSE277nnn/GSE277165/matrix/GSE277165_series_matrix.txt.gz"
gunzip "GSE277165_series_matrix.txt.gz"
