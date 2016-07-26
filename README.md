# RP-PCA
This is RP-PCA codes implemented in R. 

In source.R file, five PCA methods (CPCA, T-PCA, ROBPCA, S-ROB, RP-PCA) are coded. ROBPCA, S-ROB are implemented as PcaHubert in rrcov package and as robpca (with skew=T) in rospca package respectively. In this file, the results objects of those two functions are renamed  without any change of results values just to make them suitable our simulation machine. This file requires libraries matrixcalc, car, rotations, pragma, robust base, rrcov and rospca. rospca can be downloaded from git hub "TReynkens/rospca”, but note that to use rospca, mrfDepth package also should be installed through https://wis.kuleuven.be/stat/robust/software. Be sure that all related packages with those packages are installed your computer.

In tools.R file, data generator and working machines are coded which is essential tool of our simulation work. 

In example.R file, the first simulation of RP-PCA paper (multimodal data setting) is implemented. To run this file, source.R and tools.R files should previously be roaded.


