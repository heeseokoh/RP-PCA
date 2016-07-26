# RP-PCA
This is RP-PCA codes implemented in R

In source.R file, five PCA methods (CPCA, T-PCA, ROBPCA, S-ROB, RP-PCA) are coded. ROBPCA, S-ROB are implemented in rrcov and rospca package, but in this file, the results objects are renamed just to be suitable our simulation machine without any change of results values. This file requires libraries matrixcalc, car, rotations, pragma, robust base, rrcov and rospca. rospca can be downloaded from git hub "TReynkens/rospca”, but note that to use rospca, mrfDepth package also should be installed through https://wis.kuleuven.be/stat/robust/software.

In tools.R file, data generator and working machines are coded which is essential tool of our simulation work. 

In example.R file, the first simulation of RP-PCA paper (multimodal data setting) is implemented. To run this file, source.R and tools.R files should previously be roaded.


