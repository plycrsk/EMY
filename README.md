# EMY
EMY Thesis Scripts
These scripts were written to analyse metabolomics data, obtained in 2016/2017, from various tissues of the African
Turquoise Killifish.

The scripts are summarised as follows:

1. Metabolome_raw_data_analysis_280817 - Array manipulation of raw data, to allow for appropriate statistical tests (i.e ANOVAs) to be run. Includes False Discovery Rate control (Benjamini Hochberg). Outputs results of ANOVA for each tissue - brain, liver, heart, muscle, serum.

2. Converting_top_hits_to_PIUMet_input_280817.py - Converts output of first script into format for an online tool, PIUMET.

3. network_analysis_all_tissue_annotated_180817.py - Takes output from PIUMET (online tool), and generates visual networks.

4. Metabolome_OTU_linear_regression_290817 -Copy1.py -  Linear regression of two datasets obtained from the same fish samples.

5. meta_x_otu_regression_plots_290817-Copy1.py - Script to plot regressions and to pull/plot interesting individual datapoints from either of the original two datasets (metabolites / OTUs)


