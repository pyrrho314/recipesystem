Recipe and Primitive Tutorial
=============================

The dataset:
------------

For the sample data we are going to use US Census 2012 County Business Patterns data, specifically the "Complete Metropolitan Area File", which is 7.6MB compressed and about 63MB uncompressed, named **cbp12msa.txt**. This is a comma separated value text file. It's located at [http://www.census.gov/econ/cbp/download/](County Business Patterns), where data from 19986-2013 are available at the time of this writing.

The first issue is this file has a **.txt** extension. Data handling classes used by the recipe system are defined in configuration "kits" where all data-specific code resides. Each kit contains configuration to add particular DataHandling object to particular file extensions. We don't want to assume "*.txt" is csv, though a particular kit to set that and if placed early in the path will override the setting in other kits.



