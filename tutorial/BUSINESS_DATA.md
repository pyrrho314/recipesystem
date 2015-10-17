Recipe and Primitive Tutorial
=============================

For the sample data we are going to use US Census 2012 County Business patterns data, "Complete Metropolitcan Area File"
which is 7.6MB compressed and about 63MB uncompressed. It's a comma separated text file. The file we are looking for is at [http://www.census.gov/econ/cbp/download/](County Business Patterns). We are using the "complete Metropolitan Area File" from 2012. Data all the way back to 1986 is available, and up to 2013 at the time of this writing. 

The first issue is this file has a .txt extension. The system support data classes defined in the configuration space (inheriting from astrodata.generaldata.GeneralData to cooperate with the recipe system infrastructure, or alternately SetrefData which is a childe of GeneralData). It chooses the class based on an extention. So rename the file to csv. In an automated environment we would make a *primitive* to perform this rename, to force loading through pandas.

