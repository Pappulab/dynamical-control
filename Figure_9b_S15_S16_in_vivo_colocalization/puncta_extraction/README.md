# Overview

This repository contains the scripts used to analyze and validate the puncta centers detected as well as determine their colocalizations. This repository uses BigFISH and heavily adapts some of its example code for use on custom data.

# Files

This repository contains the following files:


```
.
├── config.ini
├── extract_and_analyze_puncta.py
├── extract_nuclei.py
├── README.md
└──  data
    ├── Nuclei_Images
    └── Colocalization_Images
├── nuclei_results
├── puncta_results
└── tools
    ├── check_summary.py
    ├── examine_titration.py
    └── export_puncta.py
```

An explanation of these files is provided here:

1. *config.ini* - This is an [INI file](https://en.wikipedia.org/wiki/INI_file) that contains the configuration for analysis. This is a companion mode of forwarding options to the main extraction script, `extract_and_analyze_puncta.py`.
2. *extract_and_analyze_puncta.py* - This is the main file that used to perform an analysis of the TIF images and extract coordinates of the puncta and their clusters. As many options are supported, many command line flags and their arguments can be passed. As that can be unwiedly, judicious defaults are set, and an option for passing in the options via an INI file is also provided.
3. *extract_nuclei.py* - This file can be used to extract the coordinates of the nuclei and their diameters.
4. *README.md* - This file.
5. *data* - This directory contains all the TIF files that will be analyzed and their data extracted.
6. *nuclei_results* - This directory contains the extracted nuclei coordinates using the `data/Nuclei_Images` dataset.
7. *puncta_results* - This directory contains the extracted puncta coordinates using the `data/Colocalization_Images` dataset.
8. *tools/check_summary.py* - This file allows the user to validate the results of the puncta and cluster extraction against that of the aggregated colocalization data.
9. *tools/examine_titration.py* - Examine the titration results and generate a summary plot.
10. *tools/export_puncta.py* - Export the puncta as txt files after the completion of *extract_and_analyze_puncta.py*.

# Execution

To execute this file with the current defaults and perform the extraction of puncta, enter this command: `python extract_and_analyze_puncta.py -d directory_containing_tif_files_here`.

Similarly, to extract the nuclei one can use this command: `python extract_nuclei.py -d data/Nuclei_Images`.