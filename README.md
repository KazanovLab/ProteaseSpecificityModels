# About

Here is method for calculating scores of proteolytic susceptibility by position-weighted matrix (PWM) of specific protease. PWMs were generated for 169 proteases from MEROPS database.

# Dependencies

For the correct work of the method, you need to install the following python libraries:
* glob
* os
* pandas
* numpy

# Quick start

1) download the archive and unzip it.
2) run the script “get_PWM_score.py”.
3) type the whole name of your input FASTA-file.
4) wait calculating PWM-scores for each proteases.
5) type the whole name of your desired output file. 

# User manual

The archive will be downloaded contains:
* The main script to get PWM scores for your sequence
* The directory “PWMs” with position-weighted matrices for 169 proteases from MEROPS database with 8 and more proteolytic events.
* The directory “annotations” with information to correspond MEROPS code to MEROPS name of a protease.
* The directory “examples” with input and output file examples.

As input you need pass classic FASTA-file where the first string is description of protein sequence, the others – sequence itself. It’s strongly recommended.

As output you will get classic CSV-file with PWM-scores for each position of protein sequence (except for last one, there you will see “Nan” value) from each protease with PWM. Information about proteases will represent as two columns: MEROPS code and MEROPS name.

# IMPORTANT

The main script, “PWMs” directory, “annotations” directory and input FASTA-file must be located in the same directory. The output file will be generated also here.

# Reporting Bugs and Feature Requests
Please use the [GitHub issue tracker](ссылка!) to report bugs or suggest features.

# Citing

