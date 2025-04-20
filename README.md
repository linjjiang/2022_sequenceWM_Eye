# Feature Misbinding Underlying Serial-order Spatial Working Memory Processing

This is an eye-tracking Study I conducted as part of my Ph.D. dissertation: [Behavioral and Neural Correlates of Serial-Order Spatial Working Memory](https://www.proquest.com/docview/3123663226?pq-origsite=gscholar&fromopenview=true&sourcetype=Dissertations%20&%20Theses)

We are preparing a manuscript to be submitted to a peer reviewed journal.

## Description

Sequential visual inputs are not processed equally. Often, we remember earlier inputs (e.g., the first item in a sequence) and more recent inputs (e.g., the last item in a sequence) more accurately than inputs in the middle, termed "primacy effect" and "recency effect" in serial-order processing literature. However, the exact neurocognitive mechanisms underlying such serial-position effects remain unknown. To investigate this, the current study used eye-tracking experiments with statistical modeling of eye-movement distributions to understand the latent neurocognitive mechanisms of serial-position effects.

In two eye-tracking experiments, we recorded eye movements from a total of 92 healthy young adults during a serial-order memory-guided saccade task. We varied how we cue the target retrieval to manipulate the involvement of different cognitive components in serial-order processing: feature misbinding vs. memory precision. A series of 2D probabilistic mixture modeling of saccade distribution were additionally conducted to analytically delineate these different cognitive components. 

Our results revealed misbinding between location and serial position as a major source of primacy and recency effects in serial-order spatial working memory as compared to varying spatial working memory precision. This suggests a redistribution of working memory resources across serial positions due to feature misbinding.

This repo contains the scripts used to perform eye-tracking preprocessing, analysis, and modeling of this project:
* preprocessing: contains scripts for preprocessing
* stats: contains scripts for aggregating preprocessed data, calculating errors/latency, checking data and modeling quality
* model: contains scripts for conducting mixture modeling
* plotting: contains scripts for replicating plots in the paper
* others: other scripts that are important for running the scripts above

## Getting Started

### Dependencies

* MATLAB R2019b or later
* [LJ_Eyes](https://github.com/linjjiang/LJ_Eyes)
* [Memfit2D Toolbox](https://github.com/johnPGrogan/MemToolbox2D)

### Installing

* Download the current repo and add it to your MATLAB PATH
```
git clone https://github.com/linjjiang/2022_sequenceWM_Eye.git
```
* Download LJ_Eyes and Memfit2D toolbox and add it to your MATLAB PATH

```
git clone https://github.com/linjjiang/LJ_Eyes.git
git clone https://github.com/johnPGrogan/MemToolbox2D.git
```

### Executing program

Under preprocessing, there are 'functions' and 'main' folder.
'main': contains main scripts for running preprocessing for each experiment and condition (exp1: long/short delay, exp2: quad/order cue)
- main_preproc_exp1_longdelay.m
- main_preproc_exp1_shortdelay.m
- main_preproc_exp2_ord.m
- main_preproc_exp2_quad.m
The following two scripts are for debugging - if you want to change the saccade selection criteria, you can change the 'select_saccades_exp1_exp2.m' 
scripts under 'functions' and then rerun the preprocessing using:
- main_preproc_exp1_longdelay_reselect.m
- main_preproc_exp1_shortdelay_reselect.m

'functions': contains important functions I wrote for preprocessing this specific task. Very very important. 
- 'cal_ord_err.m': calculate order errors
- 'check_sgolay.m': check filtering results
- 'detect_breakfix.m': detect smgs trials that break fixation during the stimulus and delay period
- 'detect_epoch_exp1_exp2.m': detect different experimental epochs
- 'load_param_exp1_exp2.m': load experimental paraemters from .dat files
- 'select_saccades_exp1_exp2.m': select and mark different types of saccades during the response window
- 'setting_exp1.m': setting file, very important. 


Under 'stats' folder, there are 'cal_stats', 'check_quality', and 'save_csv' folders
- 'cal_stats': calculate errors and other stats, including four scripts:
	- 'stats_exp1_longdelay_new.m': exp1, long delay condition
	- 'stats_exp1_shortdelay_new.m': exp1, short delay condition
	- 'stats_exp2_ord.m': exp2, order cue
	- 'stats_exp2_quad.m': exp2, quadrant cue
	- You need to run these scripts to be able to do modeling.
	
- 'check_quality': check data quality as well as goodness of fit for models:
	- 'check_quality_exp1.m': check data quality for exp1
	- 'check_quality_exp2.m': check data quality for exp2
	- You need to run the modeling before checking data quality.

- 'save_csv': save data into csv files, which can be input into other statistical softwares for further analysis
	- 'save_csv_sp_cue_exp1.m': exp1
	- 'save_csv_sp_cue_exp2.m': exp2
	- The following csv files will be saved to the output folder, e.g., './analysis/exp2'
		- p_err: primary saccade err
		- pt_err: error of primary saccade to target
		- a_err: error of closest saccade to target
		- p_rt: primary saccade latency
		- pt_rt: latency of primary saccade to target
		- a_rt: latency of closest saccade to target
		- p_orderr: order (serial position) error of primary saccade
		- p_numtr: number of trials per condition used for modeling
		- p_sd: sd of modeling output
		- p_alpha: target rate of modeling output
		- p_beta: swapping rate of modeling output
		- p_gamma: guessing rate of modeling output
		- p_aic: aic of model
		- p_bic: bic of model
		- all_err_rt: mean errors and latencies of each type of saccade of each participant

Under 'model' folder, there are two scripts for conducting modeling:
- 'model_new_exp1.m': exp1
- 'model_new_exp2.m': exp2


Under 'plotting' folder, there are four scripts for plotting figures for the paper:
- 'Fig2_exp1.m'
- 'Fig3_exp1.m'
- 'Fig4_exp2.m'
- 'Fig5_exp2.m'

Under 'others' folder, there are multiple general scripts that serves as important function for running the scripts above.



## Help & Potential Issues


## Authors

Linjing Jiang
[@linjjiang](https://github.com/linjjiang)

## Version History

* 0.1
    * Initial Release

## License

This project is licensed under the MTI License - see the LICENSE.md file for details

## Acknowledgments

* [Memfit2D Toolbox](https://github.com/johnPGrogan/MemToolbox2D)
