The following files allow one to reproduce analyses in the manuscript entitled "Alien invasion effect on soil biota: the interplay between alien species and environmental variation".

DATA & FILE OVERVIEW

***In Data folder***

The experimental data are stored in Figshare [![DOI](https://zenodo.org/badge/DOI/10.6084/m9.figshare.28139549.svg)](https://doi.org/10.6084/m9.figshare.28139549.v1).
Before the manuscript is officially published, experimental and analytical data must remain confidential. 
If needed, please contact the first or corresponding author in advance to obtain the relevant experimental data. 
All data will be made available upon acceptance of the manuscript.

*List of experimental data files (.xlsx)*

    * 1. Field_ASVs_row_data.xlsx  * 2. Field_data_group.xlsx  * 3. Field_fungi_Flattening.xlsx
    * 4. Greenhouse_ASVs_row_data.xlsx  * 5. Greenhouse_data_group.xlsx  *6. Greenhouse_fungi_Flattening.xlsx  
    * 7. Figure_S5.xlsx  * 8. Table_S6.xlsx
    
*List of phylogenetic tree data files (.newick)*  

    * 1. Field_tree.nwk  (ASVs for field survey)
    * 2. Greenhouse_tree.nwk  (ASVs for greenhouse experiment)
    * 3. IQ_tree_plant_2025.newick  (Phylogenetic tree of studying plant species)

***In Code folder***

The names of R-scripts correspond to the statistical analysis and visualization of the corresponding figures or tables in this manuscript.

*List of R-scripts*

    * 1. mian_Figure_1.R  * 2.mian_Figure_2.R  * 3. mian_Figure_3 & Table_S9.R  * 4. mian_Figure_4 & Figure_S2 & Table_S7.R  
    * 5. Table_S4 & Table_S8.R  * 6. Table_S5 & Table_S6.R  * 7. stepAICc.R (function for model selection)
    
**Data-specific onformation for:** ***Field_ASVs_row_data.xlsx***
    * Abundance table of raw sequencing data of rhizosphere fungi from field survey (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Field_data_group.xlsx***

    Variable list:
    * Sample_ID	 Sample id of  plant rhizosphere soil 
    * Latitude	 Latitude of sampling point
    * Longitude	 Longitude of sampling point
    * Year	         Sampling year
    * Site	         Name of sampling site
    * 物种名	 Chinese name of study species
    * Species        Latin name of study species
    * Genus	         Genus name of research species
    * Family	 Family name of research species
    * Origin	 Geographical origin of plants (native vs. exotic)
    * Site_pool	 The total richness of fungi for each site and year
    * Soil_ph	 Soil pH
    * Wcont	         Soil water content
    * Soil_N	 Soil total nitrogen content
    * Tave	         Mean annual temperature (℃)
    * Prec	         Mean annual precipitation (mm)
    * CV_Tave	 Coefficient of variation of mean annual temperature
    * CV_Prec	 Coefficient of variation of mean annual precipitation
    * Chol	         Leaf chlorophyll
    * SLA	         Specific leaf area (cm2 g-1)
    * LDMC	         Leaf dry matter content (g g-1)
    * SRL	         Specific root length (cm2 g-1)
    * FRR	         Fine-to-total root mass (g g-1)
    * RMF	         Root mass fraction (g g-1)
    * Funct_pPC1	 The first PCA axes of plant traits
    * Funct_pPC2	 The second PCA axes of plant traits
    * Funct_pPC3	 The third PCA axes of plant traits
    * Phylo_vPC1	 The first PCA axes of plant phylogeny
    * Phylo_vPC2	 The second PCA axes of plant phylogeny

**Data-specific onformation for:** ***Field_fungi_Flattening.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from field survey (rarefied to minimum sample size)
      
**Data-specific onformation for:** ***Greenhouse_ASVs_row_data.xlsx***
    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (not rarefied to minimum sample size)

**Data-specific onformation for:** ***Greenhouse_data_group.xlsx***

    Variable list	 Description
    * Sample_ID	 Sample id of  plant rhizosphere soil 
    * Chinese_name	 Chinese name of study species
    * Species	 Latin name of study species
    * Genus	         Genus name of research species
    * Family	 Family name of research species
    * Repeats	 Repeat number of soil samples
    * Origin	 Geographical origin of plants (native vs. exotic)
    * Hmax	         Plant height (cm)
    * Chol	         Leaf chlorophyll
    * LA	         Individual leaf area (cm2)
    * SLA	         Specific leaf area (cm2 g-1)
    * LDMC	         Leaf dry matter content (g g-1)
    * SRL	         Specific root length (cm2 g-1)
    * FRR	         Fine-to-total root mass (g g-1)
    * RMF	         Root mass fraction (g g-1)
    * Funct_pPC1	 The first PCA axes of plant traits
    * Funct_pPC2	 The second PCA axes of plant traits
    * Funct_pPC3	 The third PCA axes of plant traits
    * Phylo_vPC1	 The first PCA axes of plant phylogeny
    * Phylo_vPC2	 The second PCA axes of plant phylogeny

**Data-specific onformation for:** ***Greenhouse_fungi_Flattening.xlsx***

    * Abundance table of raw sequencing data of rhizosphere fungi from greenhouse experiment (rarefied to minimum sample size)
    
**Data-specific onformation for:** ***FungalTraits.xlsx***

    * The FungalTraits database (Põlme, S., Abarenkov, K., Henrik Nilsson, R., Lindahl, B.D., Clemmensen, K.E., Kauserud, H., et al. 2021. "FungalTraits: a user-friendly traits database of fungi and fungus-like stramenopiles." Fungal Diversity 105: 1-16.)
