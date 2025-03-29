# Coupling in situ and remote sensing data to assess α- and β-diversity over biogeographic gradients

## Description

This repository contains all the material needed to process the data, perform 
the analysis and produce the figures that can be found in 
[ [1]](https://arxiv.org/abs/2404.18485). This study compares different source of 
information (in situ and remote sensing)
to measure alpha-diversity and beta-diversity.

The repository consists of two folders: the **Data** folder, which contains 
four folders, two R scripts and one Python script 
used to store, clean and process the 'raw' data and
the **Analysis** folder containing six R scripts developed to analyze the data
and produce the various metrics, figures and tables presented in the paper. 

## Data

The **Data** folder contains four folders.

**0_Countries** contains three shapefiles used to generate the Figure 1.

**1_Plant** contains the co-occurence matrix (with sites as rows and plant 
species as columns) and the shapefile of sites (grid cells). The script 
***Extract_abc_ABC.R*** should be used to generate the file **abcABC_Plant.csv**
that will be used to compute the plant diversity metrics. 
Be patient, it can take several hours...

**2_Spectral** contains the co-occurence matrix (with sites as rows and spectral 
species as columns). The script 
***Extract_abc_ABC.R*** should be used to generate the file **abcABC_Spectral.csv**
that will be used to compute the spectral diversity metrics.

**3_Bipartite** contains the script ***Extract_Bipartite.py*** that should be
used to generate the plant-spectral bipartite network stored in the file 
**Bipartite.csv**. You can run the code using the command:

**python3 Data/3_Bipartite/Extract_Bipartite.py**

## Analysis

The **Analysis** folder contains six scripts developed to analyze the data, 
to extract the metrics and to produce the figures and tables. The scripts 
should be run successively.

***0_Map_stats_and_Filters.R*** is used to produce the Figure 1, 
Figure S1, S2, S3 and S4. It also extract the filtering information stored in
the **Filters** subfolder (created at this step).

***1_Extract_Metrics.R*** is used to extract the alpha-diversity and 
beta-diversity 
metrics stored in the **Metrics** subfolder (created at this step).

***2_Global_Analysis.R*** is used to produce the Figure 2.

***3_Global_Analysis_Filters.R*** is used to produce the Table 1.

***4_Extract_Bioregions.R*** is used to produce the Figure S5 and extract the
'bipartite communities' stored in the **Bioregions** subfolder
(created at this step). It then used to extract the bioregions (stored in
the **Bioregions** subfolder) and to produce the Figure S6, Figure 3 and Table
S1.

***5_Bioregional_Analysis.R*** is used to produce the Figure 4, Figure 5, 
Table S2 and Table S3.

## Reference and citation

If you use this code, please cite the following reference:

[1] Lenormand M, Féret JB, Papuga G, Alleaume S & Luque S (2025)
[Coupling in situ and remote sensing data to assess α- and β-diversity over biogeographic gradients.](https://nsojournals.onlinelibrary.wiley.com/doi/full/10.1111/ecog.07479) 
*Ecography* (early view).  

If you need help, find a bug, want to give me advice or feedback, please contact me!

## Repository mirrors

This repository is mirrored on both GitLab and GitHub. You can access it via the following links:

- **GitLab**: [https://gitlab.com/maximelenormand/plant-spectral-diversity](https://gitlab.com/maximelenormand/plant-spectral-diversity)  
- **GitHub**: [https://github.com/maximelenormand/Plant-spectral-diversity](https://github.com/maximelenormand/Plant-spectral-diversity)  

The repository is archived in Software Heritage:

[![SWH](https://archive.softwareheritage.org/badge/origin/https://github.com/maximelenormand/Plant-spectral-diversity/)](https://archive.softwareheritage.org/browse/origin/?origin_url=https://github.com/maximelenormand/Plant-spectral-diversity)
