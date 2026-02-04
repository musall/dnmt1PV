# DNMT1 PV data analaysis
This repo contains a Matlab-based code collection for  analysis of electrophysiological data in the study by Linde*, Yildiz*, et al 2025.
DNMT1 Coordinates PV Interneuron–Glia Coupling to Maintain Cortical Network Stability and Regulate Behavior
https://doi.org/10.1101/2025.11.18.689032

The data for the analysis can be found here: https://doi.org/10.6084/m9.figshare.31245388

The dataset contains Neuropixels recordings from 4 mice that are either PV-Cre or PV-DNMT1-KO animals. Each recording contains pre-processed LFP and spiking data that was used to create the electrophysiological results of the study.

To get to the results in the paper, download the dataset with the data in a common base folder.
You can then use the main functions (all starting with dpv_) to create different figures. For the code to work, the variable 'opts.localPath' needs to be changed to the base folder where the downloaded data was placed

Here is a short description of the included functions:


dpv_checkOscillatoryPower - Show results for analysis of LFP data as shown in Fig. 7 of the manuscript

dpv_showVisualResponse - Show current source density analysis in S1 recordings in response to visual stimulation

dpv_sensoryStimulation - Show results of spiking analysis in response to visual stimulation

dpv_pvStimulation - Show neural response to optogenetic stimulation of PV neurons
