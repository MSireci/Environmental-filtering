1) DATA: The code has been written to analysis  public available datasets, for details on them see the first section of the SI. 
They can be found in the crosssecdata.RData and longitudinal.RData files.
The file 97_otus.tree is necessary to obtain the phylogenetic distances between OTUs.

2) ANALYSIS:
The file correlation_analysis.r contains all the functions necessaries for the statistical analysis

3) EXAMPLE: to illustrate the analysis in example_analysis_glacier.r we present how to carry it out for the concrete case of the Glacier dataset (the smallest one that we have studied)

The main point is, for each biome, joining the abundances OTU datasets, to produce the pearson correlation coefficients, and the phylogenetic ones. 
The analysis of the OTU's abudances can be carried out fully in the example code presented in the folder starting by the  datasets that contains the data for many different biomes (cross-sectional or longitudinal data). 
The phylogenetic analysis has been carried out  first by producing the phylogentic trees for the OTU present in each community starting by the file  97_greengess.tree. 
Given that the production of the phylogenetic tree is a very long and tediouos work I have charged here already the output distance file for the Glacier biome to run the example of the analysis : https://www.dropbox.com/scl/fo/gzv2j64vt80p3a2afh0h7/AI-yMNHbbYUO3ttku9lflHI?rlkey=jknrigl6453vqn0wfay1fxe4y&st=qhkfh73s&dl=0

As reported in the example, the main part of the analysis consist of joining the abundances  pearson correlations with phylogenetic distances and finally obtain the average relation between the two.
In the example we also include how to randomize the phylogenetic distances to obtain a null expectations for the pattern.

