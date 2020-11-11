# TabletTVA
Encouraging digital technology in neuropsychology: The Theory of Visual Attention on tablet devices

## Overview
This repository contains the scripts and datasets needed to reproduce the TVA results in the paper "Encouraging digital technology in neuropsychology: The Theory of Visual Attention on tablet devices". 

## Contents
* Paradigm
  * CombiTVA_Asus:  Desktop paradigm
  * CombiTVA_Nexus: Tablet paradigm
* Analysis code
  * Functions: scripts for loading and pre-processing the raw data and quality checking, as well as visualisations and statistical analysis
  * 2020_ms_revision: main script creating results and visualisations for the manuscript
* Data_2020ms
  * tablettva_rawdata: raw accuracy data with added data cleaning columns 
  * tablettva_fit_fixed: LibTVA fitted data
  * tablettiva_fit_splitA and tablettiva_fit_splitB: LibTVA fitted data for split-half correlations

## Software requirements
The Matlab scripts have been tested in Matlab 2016b. The LibTVA toolbox can be downloaded from: http://www.machlea.com/mads/libtva.html
The experimental paradigms have been written and tested in Unity version 2019.1.11f1. The Tablet paradigm was performed on Android 5.1.

## Citation
For usage of the scripts and the associated manuscript, please use the following:

Tianlu Wang, Hella Thielen, Erik De Preter, Signe A. Vangkilde, Céline R. Gillebert (under review). Encouraging digital technology in neuropsychology: The Theory of Visual Attention on tablet devices.
