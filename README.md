
Authors: María Piotto, Iván Barberá, Mariano Sironi, Victoria J. Rowntree, Marcela M. Uhart, Macarena Agrelo, Alejandro A. Fernández Ajó, Jon Seger, Carina F. Marón.

Abstract
Reports of seabirds attacking marine mammals have become more frequent in the 2000s. Southern right whales (Eubalaena australis) off Península Valdés, Argentina, have been harassed by kelp gulls (Larus dominicanus) at least since the 1970s. In 2003-2013 this population experienced nine years of high calf mortality, including the highest ever recorded for a baleen whale. Using 25-year datasets (1995-2019), we studied long-term changes in gull attacks and explored whether and how they may have contributed to calf mortality. Applying generalised linear models, we found that the intensity and frequency of attacks increased sharply from 1995 to the 2000s, and that attacks on calves were three times higher than those on mothers in 2004-2019. In addition, higher intensity and frequency of gull attacks contributed to increased calf mortality, such that a year of average overall harassment would be expected to cause about twice the mortality of a hypothetical year with no attacks. Many deaths appear to occur near the end of the calving season, when older calves have suffered months of harassment. These findings suggest that chronic seabird micropredation reduces the growth of this right-whale population, and could become a major threat for other marine-mammal populations now beginning to experience seabird attacks.

This folder contains the codes to perform the analysis of the paper entitled "Seabird attacks contribute to calf mortality in a whale population".
Questions can be addressed to mpiotto@unc.edu.ar or ivanbarbera93@gmail.com

Code/Software
Here we provide R code for modelling:
(1) GAPC, GAPM and GAF over time and space (3 scripts);
•	GAPC-GAPM-GAF joint model: this is the code of the joint model. You will need the GAF data.csv and the GAP data.csv to run the model.
•	GAPC-GAPM-GAF joint model - STAN: the STAN version of the previous model.
•	Fig. 2 - Fig. S1 - Fig. S2: with this code you will be able to reproduce the Fig. 2, Fig. S1 and Fig. S2. You will not need to run the model first to reproduce these figures. However, you will need the GAF and GAP datasets used to fit the model, and also the GAF and GAP prediction table.csv and the GAP-GAF joint model_summary. For details on how these two files were created go to the script of the GAPC-GAPM-GAF joint model.

(2) Calves’ probability of dying as a function of gull-attack indexes, SST anomalies and the gulf (1 scrip);
•	Calf mortality and gull attacks: here you will find the code to run the three models on calf mortality and kelp gull attacks. You will also find how each prediction was estimated (more detail on this is available in the Extended Methods of the Supplementary Information). Codes to reproduce Fig. 3, Fig. S3 and Fig. S5 are also available in this script. You will only need the Calf mortality data.csv to run this script.

(3) The average month of calf death as a function of the three gull-attack indexes and the gulf. Also the two Pearson correlation performed (3 scripts).
•	Month of death models: we modelled the month when calves died as a function of the gulf and each of the three gull-attack indexes using this code. You will need two of our datasets to run this script: Deaths month gulf year.csv and Caf mortality data.csv.
•	Month of death plots: here you will find the code to reproduce Fig. 4A without having to run the previous script (Month of death models). To run this code you will need two files (Deaths month gulf year.csv and Caf mortality data.csv.) and also two .rsd files: month_of_death_models_prediction_table.R and month_of_death_models_residuals.rds. These .rds files are available in the repository, but you could also create both of them running the Month of death models script.
•	Dead calves' age vs. months: Pearson correlations to assess whether the month of calf death was associated with 1) the monthly average length of dead calves and 2) the relative frequency of healed and open umbilical cords. You will also find the code to reproduce Fig. 4B and Fig. 4C. The file Dead calves age.csv contains all the data that you will need to run this script.
Every script starts explaining the content of the files that you will need to run it. In addition, every script has notes and comments indicating how to use the codes and reproduce the results. However, we strongly suggest your reading the Extended Methods section of the Supplementary Information for further details on the models and to understand how some of the predictions were made.

Further questions can be addressed to mpiotto@unc.edu.ar or to ivanbarbera93@gmail.com
ome of the predictions were made.
Further questions can be addressed to mpiotto@mi.unc.edu.ar or to ivanbarbera93@gmail.com
