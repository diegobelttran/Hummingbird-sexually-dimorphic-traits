# Hummingbird-sexually-dimorphic-traits
Data and code to replicate analyses from Beltrán, Araya-Salas, Parra, Stiles, and Rico-Gevara. Manuscript submitted

This repository includes:

Code_ProcB: R code to replicate MCMCglmm models

dichromatism_all.txt: Overall dichromatism and dichromatism by patch for 237 species
Columns are -> Species:Species names,Male_chro_complexity: Interpatch male chromatic distance, Male_achro_complexity: Interpatch male achromatic distance,
Female_chro_complexity: Interpatch female chromatic distance, Female_achro_complexity: Interpatch female achromatic distance,
Dichromatism_chro: Average chromatic dichromatism, Dichromatism_achro: Average achromatic dichromatism,
Crown_dic-Gorget_dic-Belly_dic-Mantle_dic-Rump_dic: Dichromatism for each individual patch

morphologyPCA.txt: PCA scores from 14 morphological measurements for males and females of 117 species

MorphologyPCA_NoTail.txt: PCA scores from 13 morphological measurements for males and females of 117 species (same as above but without tail length)

AcousticPCs.txt: PCA scores of male song properties for 262 species

recording_metadata.csv: Metadata of recordings used in analyses

Dimorphism_Simulation.r: Script and details of dimorphism simulation


----

Tail length data from Clark (2010) can be accessed here: https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.13881

Strata from Parker III et al. (1996) can be accessed here: https://figshare.com/articles/dataset/Ecological_and_Distributional_Databases_for_Neotropical_Birds/8956955

Min, mid, and max elevation data from Rangel et al. (2015) can be accessed here: https://onlinelibrary.wiley.com/doi/abs/10.1111/evo.12644
----

McGuire et al. (2014) MCC tree (293 spp.), 
Random sample of 100 trees from posterior distribution from McGuire et al. (2014),
Weight-corrected (by square root) tail length dimorphism for 172 species using Clark (2010),
and specific strata and altitude hummingbird data used in this study will be happily shared upon request (diegofbelttran@gmail.com)