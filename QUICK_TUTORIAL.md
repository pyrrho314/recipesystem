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

    $ pwd
    /home/craiga/nrm_demo
    $ git clone https://github.com/pyrrho314/nrm_sample_data.git
    Cloning into 'nrm_sample_data'...
    remote: Counting objects: 9, done.
    remote: Compressing objects: 100% (8/8), done.
    remote: Total 9 (delta 0), reused 6 (delta 0), pack-reused 0
    Unpacking objects: 100% (9/9), done.
    Checking connectivity... done.

Make a working directory
-------------------------

Recipe processes generally output data in the current working directory, but can load data from some other input directory.  This ensures input data is not overwritten. Further, in this model, data objects descended from the SetrefData class (in kit_Novem) will always write data into the current directory by default. I.e. if you load a data sets at "/data/raw_data/the_data.xls", and then immediately writes it without setting the filename, it will be written to ./the_data.xls. This makes it easy to to delete output and begin again without harming input sources.

In ~/nrm_demo make a working directory:

    $ mkdir tutorial_wd
    $ cd tutorial_wd
    
To execute a simple recipe:

    $ kit -r showInputs ../nrm_sample_data/*.xls
    
This should output the following (albeit in color, to turn off color either use the kit command line option (see`kit --help`) or export NOCOLOR=nocolor:

    $ kit -r showInputs ../nrm_sample_data/*.xls
    About to process 1 lists of datasets.
    Starting Reduction on set #1 of 1
      Processing single dataset:
        /home/craiga/nrm_demo/nrm_sample_data/6-digit_2012_Codes.xls  
    ================================================================================
    PRIMITIVE: showInputs
    ================================================================================
    1 inputs
    #1
    filename  : /home/craiga/nrm_demo/nrm_sample_data/6-digit_2012_Codes.xls
    data types: ['SETREF', 'TABLE', 'UNINGESTED']
    data_obj  : <class 'pandasdata.PandasData'>
    (r1048) End of Recipe: WRITING FINAL OUTPUTS:-----------------------------------
    Writing from Context Stream (r1048): main : 1 datasets in stream                
    file: /home/craiga/nrm_demo/tutorial_wd/6-digit_2012_Codes.xls
          Does not need writing.
    (r1111) End of Recipe: WROTE FINAL 0 OUTPUTS of 1                               

If this works the system is operational,  See [[QUICK_TUTORIAL_BUSINESS_DATA.md]].

