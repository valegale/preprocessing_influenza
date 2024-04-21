# Preprocessing Scripts for Influenza virus files

This repository contains a preprocessing pipeline for Influenza proteins data file. The pipeline offers several options for preprocessing data, including removing reassortant entries, translation, and removal of duplicates.

### The data should be organized as follows:

{path_folder}/

│

├── ha_cleaned/
│ ├── ha_{bin}.fa
│ └── ...
│
└── na_cleaned/
├── na_{bin}.fa
└── ...

markdown
Copy code

- `{path_folder}`: This is the main folder containing the data.
- `ha_cleaned/`: This subfolder contains the hemagglutinin (HA) sequences after initial manual curation.
  - `ha_{bin}.fa`: This file contains HA sequences for a specific flu season bin (e.g., `ha_10_11.fa` for the flu season 2010-2011).
- `na_cleaned/`: This subfolder contains the neuraminidase (NA) sequences after initial manual curation.
  - `na_{bin}.fa`: This file contains NA sequences for a specific flu season bin (e.g., `na_10_11.fa` for the flu season 2010-2011).

The sequences in these files have undergone initial cleaning steps, such as removing partially sequenced strains and sequences with a large number of ambiguous characters, which are not included in these preprocessing scripts.

Additionally, two more variables can be modified at the beginning of the file:

- `first_bin`: The first flu season bin to analyze.
- `last_bin`: The last flu season bin to analyze.

These variables allow you to specify a subset of files to process if needed. For examp


Below are the available options:

## Options

### 1. Remove Reassortant

This option removes reassortant entries from the data files. Reassortants are entries that contain mixed genetic material and are often filtered out in certain analyses.

### 2. Translate

This option translates the data files into another language. Translation can be useful for multilingual analysis or for adapting the data to a specific audience.


### 3. Remove Duplicates

This option removes duplicate entries from the data files. Duplicate entries can skew analysis results and are often removed to ensure data integrity.

### 4. Full Pipeline (Not Suggested)

This option runs the full preprocessing pipeline, including removing reassortants, translation, and removing duplicates. However, using this option is not suggested as it is better to run the three steps separatedly and check the intermediate results before proceeding to the next step (first remove the reassortant, then translate, and finally removing duplicates).
