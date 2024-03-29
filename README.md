This repository contains R and MATLAB Code for our paper "EEG-based clusters differentiate psychological distress, sleep quality and cognitive function in adolescents"

Full paper here: http://bit.ly/EEG-clusters

Summary blog here: https://www.linkedin.com/pulse/can-simple-test-measures-brain-activity-help-identify-owen-forbes/?trackingId=0Ax6NF7vSi2qv%2FdjHE9Vkg%3D%3D

* preprocessing_pipeline_matlab_public.m is MATLAB code for an automated pre-processing pipeline for cleaning and artifact removal etc from resting state EEG data

* multitaper_matlab_public.m is MATLAB code for multitaper analysis of cleaned EEG data

* EEG_clusters_public.Rmd is R code for our multi-stage analysis pipeline to calculate EEG frequency features, do dimensionality reduction, apply multiple clustering methods, and test for differences in health & cognitive function measures between EEG-based clusters
