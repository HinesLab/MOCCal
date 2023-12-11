# MOCCal

MOCCal, or Multi-Omic CCS Calibrator, is a Python application for high-dimensional, multi-omic traveling-wave ion mobility mass spectrometry (TWIM-MS) data anlaysis. Functionality includes collision cross section (CCS) calibration, experimental data biomolecular class assignment, and experimental class-specific CCS calculations. Notably, MOCCal offers class assignment and CCS calculations without need for identifying the features first. 

For installation instructions, please see below.
For python scripts, usage tutorials, data templates, and example data, please see the UserData or RawData folders.

## Instructions for running MOCCal (UserData)

If you plan to use MOCCal with processed calibration data, you do not need to install Python or any of the packages required to run the code. To use this version of MOCCal, simply download the MOCCal zip file (https://outlookuga-my.sharepoint.com/:f:/g/personal/kmh84885_uga_edu/Evrw4fZFIRNAsAEXZu7DqrABvW48aiGqHaXCUUVFy3iD8g?e=BUxa4s), extract all, and use the MOCCal executable file. The only requirement is that the MOCCal executable is in a location that also contains a folder titled “Output.” The python script from which the executable was built and a guide and data templates/examples are in the UserDT folder.

## Instructions for running MOCCal (RawData)

If you plan to use raw calibration files, you will need to install pnnl’s DEIMoS, version 1.3.2 (http://github.com/pnnl/deimos). After DEIMoS is set up, you can then run MOCCal_RawData.py in the DEIMoS virtual environment. Data templates/examples are in the RawData folder. Raw data example files can be found in the RawData_ExampleData zip file at https://outlookuga-my.sharepoint.com/:f:/g/personal/kmh84885_uga_edu/Evrw4fZFIRNAsAEXZu7DqrABvW48aiGqHaXCUUVFy3iD8g?e=BUxa4s

## Citing MOCCal

If you would like to reference MOCCal in an academic paper, we ask that you include the following:

•	MOCCal https://github.com/HinesLab/MOCCal (accessed MM/YYYY)

•	PUBLICATION CITATION

## Getting Help

If you have any questions about this software, please email Kelly Hines (Kelly.Hines@uga.edu) or post to our GitHub page (https://github.com/HinesLab/MOCCal) .  

