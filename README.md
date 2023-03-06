# Death by a Thousand Cuts – Combining Kinase Inhibitors for Selective Target Inhibition and Rational Polypharmacology
Ian R. Outhwaite, Sukrit Singh, Benedict-Tilman Berger, Stefan Knapp,  John D. Chodera, Markus A. Seeliger\
doi: https://doi.org/10.1101/2023.01.13.523972

## Running the MMS Method
Download the files in run_MMS/ to a working directory.\
Instructions for running the MMS method are in the included in the `run_MMS/readme.md` file.\
Please report any issues or problems here or to ian.outhwaite@stonybrookmedicine.edu

## Installation
Install the necessary dependencies using the Anaconda or Mamba package manager\
To install run: `conda env create -n mms --file environment.yaml`

**Abstract** \
Kinase inhibitors are successful therapeutics in the treatment of cancers and autoimmune diseases and are useful tools in biomedical research. The high sequence and structural conservation of the catalytic kinase domain complicates the development of specific kinase inhibitors. As a consequence, most kinase inhibitors also inhibit off-target kinases which complicates the interpretation of phenotypic responses. Additionally, inhibition of off-targets may cause toxicity in patients. Therefore, highly selective kinase inhibition is a major goal in both biomedical research and clinical practice. Currently, efforts to improve selective kinase inhibition are dominated by the development of new kinase inhibitors. Here, we present an alternative solution to this problem by combining inhibitors with divergent off-target activities. We have developed a multicompound-multitarget scoring (MMS) method framework that combines inhibitors to maximize target inhibition and to minimize off-target inhibition. Additionally, this framework enables rational polypharmacology by allowing optimization of inhibitor combinations against multiple selected on-targets and off-targets. Using MMS with previously published chemogenomic kinase inhibitor datasets we determine inhibitor combinations that achieve potent activity against a target kinase and that are more selective than the most selective single inhibitor against that target. We validate the calculated effect and selectivity of a combination of inhibitors using the in cellulo NanoBRET assay. The MMS framework is generalizable to other pharmacological targets where compound specificity is a challenge and diverse compound libraries are available.

**Competing Interest Statement** \
B.-T.B. is the CEO and a shareholder of CELLinib GmbH, Frankfurt, Germany. J.D.C. is a current member of the Scientific Advisory Boards of OpenEye Scientific Software, Interline Therapeutics, and Redesign Science. The Chodera laboratory receives or has received funding from the National Institute of Health, the National Science Foundation, the Parker Institute for Cancer Immunotherapy, Relay Therapeutics, Entasis Therapeutics, Silicon Therapeutics, EMD Serono (Merck KGaA), AstraZeneca, Vir Biotechnology, XtalPi, Interline Therapeutics, and the Molecular Sciences Software Institute, the Starr Cancer Consortium, the Open Force Field Consortium, Cycle for Survival, a Louis V. Gerstner Young Investigator Award, and the Sloan Kettering Institute. A complete funding history for the Chodera lab can be found at http://choderalab.org/funding. No other authors declare competing interests.

**This work is currently in the preprint stage and has not yet been peer reviewed.**
