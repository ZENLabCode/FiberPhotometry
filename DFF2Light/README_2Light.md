# Fiber Photometry GCaMP Data Analysis

This repository contains scripts for the analysis of Fiber Photometry GCaMP data before, during, and after light pulse periods.

If you have any questions or require the relevant publication to cite for this code, please contact:
- Micaela Borsa: micaela.borsa@unibe.ch
- Antoine Adamantidis: antoine.adamantidis@unibe.ch
	
## Workflow Overview

1. Start with `Cut_lvm2light_DFF.m` to analyze all individual raw files:
	- This script calculates the ΔF/F (delta-fluorescence over baseline-fluorescence) of your photometry signal per experiment.
	- The ΔF/F is cropped according to the "lights-on / lights-off" periods, which are recorded in a reference channel (refCHA).

2. Prepare the `filelist.xlsx` file:
	- Create a file listing the paths of the Ch*.mat files along with the following information: Group; Trial; Condition; Channel (see template in the folder).

3. Run `Analyze_DFF_Grouped_2Light.m`:
	- This script analyzes the mean ΔF/F data for each animal and condition for specific time profiles: before, during, and after the light pulse.



## Notes

The scripts are currently tailored for internal and specific use; they are not fully automated or optimized for a broader audience.
Some adaptability is however possible in the "Parameters" section (first one) of each script.

## Pubblication
These scripts have been used for the following pubblication:

Brenna Andrea, **Borsa Micaela**, Saro Gabriella, Ripperger Jürgen A, Glauser Dominique A, Yang Zhihong, **Adamantidis Antoine**, Albrecht Urs (2024) *Cyclin-dependent kinase 5 (Cdk5) activity is modulated by light and gates rapid phase shifts of the circadian clock*, eLife 13:RP97029; https://doi.org/10.7554/eLife.97029.2


