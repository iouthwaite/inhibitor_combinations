Analysis:

1. Removed all kinases from the three datasets (Fabian et al 2005, Karaman et al 2008, PKIS2) except the EPH family kinases

2. Ran JS_singletarget.py using the following settings for each dataset:

Prior: Poisson with either mu=500 or mu=1000 in order to better capture lower off-target effects in the data
Prior size: 10,000
noise variance: 2.5
max combination iter: 5 if Fabian or Karaman, 3 if PKIS2
influence: 0

Set the minimum inhibitor potency to 70 in order to facilitate analysis

The step for weight optimization of inhibitors was a factor of: 4

Lines numbers of changes:
131, 167, 302

Conducted five technical replicates for each setting

Used the script "run.py" in order to produce the results
