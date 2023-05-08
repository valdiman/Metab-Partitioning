# Partitioning

<a rel="license" href="http://creativecommons.org/licenses/by/4.0/"><img alt="Creative Commons License" style="border-width:0" src="https://i.creativecommons.org/l/by/4.0/88x31.png" /></a><br />This work is licensed under a <a rel="license" href="http://creativecommons.org/licenses/by/4.0/">Creative Commons Attribution 4.0 International License</a>.

----------------------
General Information
----------------------

Deposit Title: PCB partitioning in multiwell setup experiments

Contributor information:

Andres Martinez, PhD
University of Iowa - Department of Civil & Environmental Engineering
Iowa Superfund Research Program (ISRP)
andres-martinez@uiowa.edu
ORCID: 0000-0002-0572-1494

This README file was generated on May 8 2023 by Andres Martinez.

This work was supported by the National Institutes of Environmental Health Sciences (NIEHS) grant #P42ES013661.  The funding sponsor did not have any role in study design; in collection, analysis, and/or interpretation of data; in creation of the dataset; and/or in the decision to submit this data for publication or deposit it in a repository.


This README file describes a revised version of scripts for calculating the individual fraction of PCBs in multiwell experiments used in this paper:

The script is developed for individual PCBs, OH-PCBs and sulfate-PCBs from experiments from Porject 1, ISRP.

--------
FILE OVERVIEW
--------

File Name: MetabPartitioning
File Size: 9 kb
File Format: .R
File description: Contains the raw code for the calculation of individual PCB fractions in
the multiwell expereriments in "R".

--------
PREREQUISITES & DEPENDENCIES
--------

This section of the ReadMe file lists the necessary software required to run codes in "R".

Software:
- Any web browser (e.g., Google Chrome, Microsoft Edge, Mozilla Firefox, etc.)
- R-studio for easily viewing, editing, and executing "R" code as a regular "R script" file:
https://www.rstudio.com/products/rstudio/download/

--------
SOFTWARE INSTALLATION
--------

This section of the ReadMe file provides short instructions on how to download and install "R Studio".  "R Studio" is an open source (no product license required) integrated development environment (IDE) for "R" and completely free to use.  To install "R Studio" follow the instructions below:

1. Visit the following web address: https://www.rstudio.com/products/rstudio/download/
2. Click the "download" button beneath RStudio Desktop
3. Click the button beneath "Download RStudio Desktop".  This will download the correct installation file based on the operating system detected.
4. Run the installation file and follow on-screen instructions. 

--------
HOW TO RUN THE CODE
--------

It is recommended to create a project in RStudio and name it Metab-Partitioning. Paste and execute Subfolders.R code to generate the different folders. Then paste the MetabPartitioning.R into the R folder just generated. All the packages needed, with their respecting libraries, are included in the code. All the needed data are in the code. If needed, the fraction results can be exported as a csv file to the Output/Data/csv folder.


