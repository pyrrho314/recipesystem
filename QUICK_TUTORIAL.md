Preparation
============

Install Recipe System and Novem Kit
------------------------------------

For installation see the [README.md](./README.md). Installation includes two repositories, 
one is the infrastructure, the `astrodata` module. The name stems from the systems initial development for
astrophysical data. However, data transformation does not happen in this module, only data handling and process 
execution.

Get the sample data
--------------------

Pull the data to the same directory used in install, for convienience, and so the relative paths in the examples work.  Otherwise adjust the pathnames accordingly.

    craiga@pleiades:~/nrm_demo$ pwd
    /home/craiga/nrm_demo
    craiga@pleiades:~/nrm_demo$ git clone https://github.com/pyrrho314/nrm_sample_data.git
    Cloning into 'nrm_sample_data'...
    remote: Counting objects: 9, done.
    remote: Compressing objects: 100% (8/8), done.
    remote: Total 9 (delta 0), reused 6 (delta 0), pack-reused 0
    Unpacking objects: 100% (9/9), done.
    Checking connectivity... done.

Make a working directory
-------------------------

Recipe processes generally output data in the current working directory, but can load data from some other input directory.  This ensures input data is not overwritten. Further, in this model, data objects descended from the SetrefData class (in kit_Novem) will always write data into the current directory by default. I.e. if you load a data sets at "/data/raw_data/the_data.xls", and then immediately writes it without setting the filename, it will be written to ./the_data.xls.

So in ~/nrm_demo make a working directory:

    mkdir tutorital_wd
    
