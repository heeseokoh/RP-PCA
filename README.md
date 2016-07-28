# RP-PCA

Here are codes for reproducing the first1 simulation introduced in RP-PCA paper.

In source.R file, five PCA methods (CPCA, T-PCA, ROBPCA, S-ROB, RP-PCA) are coded. In fact, functions of ROBPCA and S-ROB are already implemented as PcaHubert in rrcov package and as robpca (with skew=T) in rospca package respectively. In this file, the resultant objects of those two functions are renamed without any change of the resultant values just to make them suitable for our simulation machine.

In tools.R file, a data generator and a working machine are coded and they are essential tools of our simulation work. 

In example.R file, the first simulation of RP-PCA paper (multimodal data setting) is implemented. To run this file, source.R and tools.R files should previously be loaded.

 These files require libraries matrixcalc, car, rotations, pragma, mvtnorm, fBasics, robustbase, rrcov and rospca. rospca can be downloaded from https://github.com/TReynkens/rospca. Note that to use rospca, mrfDepth package also should additionally be installed through https://wis.kuleuven.be/stat/robust/software. Be sure that all other related packages with those packages are installed on your computer before you run the files.


