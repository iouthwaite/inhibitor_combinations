Running the MMS Method

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

0. Download dependencies

Use the enviroment.yaml located in /inhibitor_combinations

These scripts are supported for python -v 3.9.13 with the following dependencies:

numpy, pandas, openpyxl, xlsxwriter, itertools, random, statistics, math, scipy, multiprocessing

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

1. Define your input data

Users can select from some preformatted datasets from prior studies:

Datasets included in /datasets and preformatted:
    
    	A. Inhibitors considered for inclusion in PKIS2, DOI: 10.1371/journal.pone.0181585 (Supplemental Table 4)
    	B. Davis 2011, DOI: 10.1038/nbt.1990 (Supplemental Table 4, human kinases only)
    	C. Fabian 2005, DOI: 10.1038/nbt1068 (Supplemental Table 4)
    	D. Karaman 2008, DOI: 10.1038/nbt1358 (Supplemental Table 2)
    	E. Klaeger 2017, DOI: 10.1126/science.aan4368 (Supplemental Table 3: Kinobeads Drugmatrix by Protein)

If a user wants to include their own data, please be sure to format it in the same way as the included datasets and input values on an activity scale

Specifically,

	A. Label columns by targets and rows by compounds
	B. Label the column for compounds with 'Compound'
	C. Leave dummy columns in the same locations as in the preformatted datasets
	D. Convert inpout data to the activity scale, given a reference concentration frame. The program output concentrations will be relative to this input concentration frame.

    	For example, given Kd values in nM and a desired reference concentration range of 1 uM
    	activity = 100% / [(Kd/1,000 (M))+1].

2. Set the parameters that you want to use

Program parameters are as follows, and can be set in the script "MMS_run.py"

	A. single, a boolean for single-target analysis (tests all targets in a dataset one by one) or multiple-target analysis (tests combinations for co-inhibiting multiple targets)

	B. dataset_path, a string defining the path to the input dataset

	C. num_inhibitors, an int denoting the number of inhibitors (from the top of the input dataset going down by column) to use in the analysis. Set to the total number to use or, or a subset if desired.

	D. target_kinases, a list of strings with the names of the targets for multiple-target analysis. Ignored for single-target analysis.

	E. non_off_target_kinases, a list of strings with the names of targets in the input dataset to NOT include when determining the off-target effects of compounds. For example, if an input dataset has data for multiple mutants, those might be considered targets, but not off-targets, since they would not normally co-occur.

	F. on_target_activity_threshold, an int from 1->99 that sets the threshold for minimum on-target activity, given the activity of compounds used in the input reference frame. For example, a threshold of 90 would mean that the program would only select compounds that reach 90% target activity when determining which inhibitors to use to generate combinations. If no compounds reach this threshold then no results will be generated for the target(s).

	G. max_number_inhibs_per_combo, an int that defines the maximum number of inhibitors to try in a combination (ex: up to 4 inhibitors in a single combination)

	H. num_UW_to_test, an int that describes the number of combinations to consider for concentration optimization and then scoring

	I. num_p, an int for the number of procs to call. Optimizing this variable and runtime will depend on the user's workstation setup.

	J. n_replicates, an int for the number of technical replicates of the analysis to perform. 

	K. outfile_prefix, a string for the prefix to use for output files

	L. lower_limit, an int for the lower limit for compound dilution(s) given input reference frame. Ex: given a reference frame of 1uM, 1000 would mean using compounds at a lower limit of 1 nM

	M. upper_limit, an int for the upper limit of the concentration of compound(s)s given input reference frame. Ex: given a reference frame of 1uM, 100 would mean using compounds at an upper limit of 100 uM

	N. max_num_to_report, an int for the number of results to report per target

ADDITIONAL SETTINGS: users will generally not change these settings.

	O. prior_type, an int for whether to use a penalty prior with a base shape of a poisson or beta distribution? Poisson is suggested

	P. bdist_alpha, set the parameters for the beta penalty prior distribution

	Q. bdist_beta, set the parameters for the beta penalty prior distribution

	R. mu_range, set the parameters for the Poisson penalty prior distribution

	S. influence, set the high off-target penalty (0.1->0.3 is suggested)

	T. noise_variance, set the noise added to off-target measurments when building off-target distributions. 2.5 is suggested, set to 0 if you want to improve reproducibility of results. This value corresponds to the variance of the normal distributions that are used for sampling each off-target calculation

	U / V : Step sizes
	set the step sizes for calculating optimal concentrations of inhibitors
	The R1 step refers to how much inhibitor pairs are increased or decreased by (a factor of R1) during initial branching
	The R2 step refers to the factor that concentrations are decreased by at each step in order to obtain minimum on-target activity
	Increasing the step sizes will increase the speed of the optimization steps, but will reduce the accuracy of obtained concentrations
	We reccomend using larger step sizes for initial screens (ex: R1=5, R2=1.2), and then decreasing the step sizes for smaller analysis using drug-target pairs
	U. R1_step
	V. R2_step

	W. time_lim, an int for setting the time limit per single target-combination in the optimization concentrations protocol (in seconds)

	X. size, an int for setting the number of points samples from the prior distribution used to build the penalty prior

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

3. Executing the program

In your terminal, execute MMS_run.py with the command "python MMS_run.py"

   	A. Be sure that the path given to the input dataset matches what is in the MMS_run.py script
   	B. Be sure that the script "MMS_run.py" is in your current directory (or that you call the relative path)
   	C. Be sure that the script "MMS_methods.py" is in the same directory as "MMS_run.py"

This will inititate the program, in that terminal. You can set the number of procs to call by adding that number after the command if desired, this will overwrite the setting in the MMS_run.py file. Ex: "python MMS_run.py 20"

The program will create several output .xlsx files that can be converted to .csv or analyzed as desired by the user.

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

4. Troubleshooting and common issues

My output files are blank for certain combination numbers. Why is this?
>There were no possible combinations of inhibitors that could be used. For example, this is common at the i=1 condition for multiple-target analysis, since you would need to find a single inhibitor that is potent for all of your targets.

My JSD scores are decent but there are high off-target effects above 95%? What is going on?
>In certain rare cases early in testing we observed that the JSD score gets stuck if the concentration steps are too high and the off-target penalty for high off-target effects is relatively low. This leads to a subset of off-targets being potently inhibited in a single bin in the off-target probability distribution which is favorable when computing the JSD metric. We've since implemented changes to how we obtimize concentrations in order to prevent this from happening, and do not observe this behavior in out implementations, but we cannot account for all user cases. If this behavior is observed, try decreasing the "R1_step" parameter size, or increasing the "influence" parameter step size to ~0.3.

I'm getting JSD scores of 1! Are my inhibitors perfectly selective?
>No this does not mean that an inhibitor is perfectly selective -- it just means that, given the range of off-targets being scored using the particular penalty prior for that analysis, there were no off-target effects in the activity range following dilution of the inhibitor(s) to minimal on-target activity. If you know that there are off-targets and want to score these with better precision, try increasing how broad the penalty prior is, for example, using a mu=1200 (or even 1700) value to better score low off-target effects.

My JSD score for a combination is higher than for a single inhibitor. Does this mean that the combination is better?
>Not necessarily. We suggest conducting multiple technical replicates in order to statistically evaluate improvements in JSD score. We also suggest using an absolute-value cutoff since very small improvements may be non-meaningful in practice. The value that a user selected will depend on their use case and the data being studied.

What is a "good" JSD score?
>This is relative, since it depends on the type of penalty prior being used and the number of possible off-targets. In general, using a more broad penalty prior (ex: mu=700 versus mu=200) will decrease the JSD score, since a greater range of off-target activities will be penalized. While it depends on numerous factors and hard cutoffs are not ideal, JSD scores above 0.9 generally reflect nicely selective compounds.

The iterations during runtime sometimes drop and then go back up. What is up with that?
>The print call is specific to the current process, which doesn't always return in the order it was generated. In general however, the iters will increase, although sometimes this number may bounce around a little, and this does not affect program function.

The script is taking a prohibitively long time to run! Help!
>The following can improve runtime
     A. Increase the R1_step size and use this first run as a screen to hone down in on combinations that might be beneficial which can be analyzed in a second smaller run with just select inhibitors and targets
     B. Incrase the number of procs called (this will not always help and depends on your workstation setup)
     C. If certain targets are very promiscuous they will have lots of possible combinations. If you don't care about these targets, removing them from the analysis will speed up processing for all the other targets.
     D. Increasing the on-target threshold will decrease the number of inhibitors considered and improve runtime, although this may not be desireable depending upon the use case
     E. Generating high combinations of inhibitors can take a long time, and depending upon the data, may not yield an improvement in the score. We suggest using n+2 over the number of targets (ex: up to three inhibitors for a single target, five for triplets of targets) at most before exploring higher combination numbers.




                   
