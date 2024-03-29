# Releases of domaps and associated webapps

## Version 1.0.1 - development in progress

### Features
 - Neighborhood analysis from EVprofiler (Martin-Jaular, et al. EMBO Journal 2021) added.

## Version 0.2 - development in progress (at revised version submission)

### Features
 - Now highly customizable data upload with gui elements as part of the new subpackage 'gui'
 - Included install and testing action on github for python 3.7-3.10.
 - PCA visualization module including elbow and loading plots
 - Compatibility with custom organelle annotation and complex member annotations
 - Expanded help texts significantly
 - SVM misclassification upload now also for MetaMass and custom input
 - Correlation heatmaps and pairwise scatter plots included

### Changes
 - Bumped version of underlying packages
 - SVM interface now more flexible and including class size as parameter
 - Intramap scatter now always normalized and different representations available

### Fixes
 - Tabs in the interface are no longer scrolling unnecessarily

### Known bugs
 - The reannotation of gene names from uniprot is currently disabled, as it needs to be adapted to the new Uniprot REST API

## Version 0.1 (at initial journal submission)

### Features

### Changes
 - SMV table upload now inside the benchmarking tab
 - Different display of intra-map scatter

### Fixes

## Version 0.1 (at biorchive submission)

### Features
 - Upload of files from MaxQuant, Spectronaut and custom sources
 - Configuration of standard data processing pipeline with parameters
 - Experiment quality control and benchmarking options
 - Separate data management section for joining, renaming and amending data collections
