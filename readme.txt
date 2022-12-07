
Automated and standardized data analysis for the Plasmodium falciparum in vitro Parasite Reduction Ratio (PRR) assay.
October 2022, Version 8.3.1

	Purpose:
	This code is designed for analysis of PRR assay data, i.e. the number of parasites surviving drug treatment (viable parasites) over time of drug exposure.
	It provides you with the lag phase of the drug, the PRR (log10 drop of viable parasites within 48 h) and the 99.9% parasite clearance time (time to kill 99.9% of the initial 100 000 parasites).

	How to use this code:
	- The input file containing the raw data (log10 viable parasites, normalized) should have the same format as provided in the example (input_example.xslx).
	- The folder containing the Code (the .R file) should further contain the following sub-folders from which data can be drawn / to which results can be pasted
	  (with the exact same spelling): 
		a. "data": this folder should contain the Excel file with the input data.
		b. "results": this folder will be filled with an Excel file containing the relevant pharmacodynamic parameters.
		c. "figures": this folder will be filled with figures (one per compound).
	- Note that the main folder must always contain the file with the R-code unless you adjust the working directory manually.
	- Before running the code, you need to specify the name of the input file in the R code in line 26.
	- When first using the code, you might need to install some of the required packages. For this, simply remove the "#" at the beginning of line 12, leave the cursor in that line and press "Run". As soon as all the packages were installed, you can re-add the "#".
	- Then you can run the code.
	
	Good to know:
	- You can run a single campaign or several campaigns at once. Just make sure the campaign identifier in the input file is not the same.
	- You must always include data for the GROWTH CONTROL (untreated samples at 0 and 48 h). If you do not have these data, you can insert mock data.

Copyright 2022 PRR assay V2 core team

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.
You should have received a copy of the GNU General Public License along with this program. If not, see <https://www.gnu.org/licenses/>.
