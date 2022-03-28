Analysis:

1. Removed all kinases from the three datasets (Fabian et al 2005, Karaman et al 2008, PKIS2) except the EPH family kinases

Ended up just using results for the PKIS2 dataset

Edited JS single target so that no inhibitor dilutions were performed (same as for PKIS2 general screen - NOT like Klaeger & Karaman)

2. Ran JS_singletarget.py using the following settings for each dataset:

Prior: Poisson with either mu=200 or mu=1200 in order to better capture lower off-target effects in the data
Prior size: 10,000
noise variance: 2.5
max combination iter: 3
influence: 0.1
Opt step sizes: 2, 1.1

Set the minimum inhibitor potency to 90

Conducted five technical replicates for each setting

Used the script "run.py" in order to produce the results
