# EyeFunctions-Python-Module
eyeFunctions is a custom-made Python module used for processing and analysing time series data. It was initially built for eye-tracking data, but was then adapted to handle python lists of up to 3 dimensions.
eyeFunctions is composed by three subpackages:
- the 'toolBox.py' module contains a rich set of functions used to handle python lists
- the 'preprocessing.py' module contains a set of functions that are necessary (and in most cases sufficient) to preprocess eye data converted in .asc, .csv, or .txt files. It can also be useful to preprocess other time series data, such as EEG data.
- the 'analyse.py' module contains a set of functions allowing to perform simple analyses for time series (repeated t-tests, permutations, etc).

If you struggle using the module or need additional information, please send your questions to sylvain.gerin@uclouvain.be
Please note that most of the functions created here are inspired from this paper: Mathôt, S., Vilotijević, A. Methods in cognitive pupillometry: Design, preprocessing, and statistical analysis. Behav Res 55, 3055–3077 (2023). https://doi.org/10.3758/s13428-022-01957-7
