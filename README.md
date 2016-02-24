recipesystem
============

The code in the *recipesystem* repository is released under the Mozilla Public License v. 2 and is meant to be used with the Novem Recipe Kit, at [https://github.com/pyrrho314/recipe_kits](https://github.com/pyrrho314/recipe_kits). Collectively these are the Novem Recipe Machine. 

The recipe system executes data transformations in a controlled way. The scientific programming environment is simple, the author does not need to understand any particular object models of the recipe system, and writes "primitive" data transformations in the form of python generator functions which, although of course any python sophistication is possible, can be written much like scripts with "yield" statements. This makes existing python code easy to translate into the recipe system, which means existing code can be automated, monitored, and controlled with minimal refactoring at first. The recipe system can also be used interactively or embedded in other projects. The code behind the `kit` command which runs recipes simply calls a factory for an object which can execute recipes. Any python program can thus use this to execute recipies and primitives. 

Having said that, a proper recipe system will not limit itself to something similar to a well controlled family of related scripts, but will have good object modeling for it's data. These objects will become sophisticated helpers handling the specific type of data involved. However, they initially can and will start as very simple, perhaps empty, classes descending from `astrodata.generaldata.GeneralData`.  This provides the minimal interface needed by the recipe system which will probably include writing a "load_header", "load" and "do_write" functions. Generally one will want also to add the get and set property functions, to help with this one can descend from the `SetRefData` family of classes in the Novem Kit configuration package.

*See installation instructions below, also a simple tutorial (sample code is included in the Novem Recipe Kit).*

brief history
-------------
The Novem Recipe Machine software is derived from the Gemini Observatory Python Package started circa 2005. The system model involves a core infrastructure (the Recipe System) which loads external "recipe kits" which contain all the code specific to a particular organization or data type. Gemini released this system at ADASS 2013 in Kona, Hawaii, under a BSD license, though currently there is no public repository for that version of the system. It is planned to arrive on GitHub or some other public git repository, and also to merge changes made to the NRM fork.

The Novem Recipe Machine version is released under MPL 2. 

At the present time, Novem LLC is generalizing the system for commercial use and using the system to automate its clients' particular data processing pipelines using the recipe model. Key development goals for the immediate future:

* relocation of AstroData Class source code into the astrodata_FITS configuration package
* general purpose documentation (available documentation is Astronomy Oriented)
* automatic cloud-based scaling, job cooperation
* general purpose HTML5/Javascript interfaces for system control and monitoring
* generic library of numpy and pandas primitives
* expansion of the general purpose Set Reference Data Type
* creation of a project roadmap

**Install**
-----------

I'm creating a directory `nrm_demo` and executing the following there.

First get the instrastructure in the Novem Recipe System from git.

    git clone https://github.com/pyrrho314/recipesystem.git

We'll add a script directory to the PATH environment variable below, and a python directory to PYTHONPATH, but first, get the standard recipe kit, kit_Novem, in the recipe_kits directory.

    git clone https://github.com/pyrrho314/recipe_kits.git
    
This gives us:
    craiga@pleiades:~/nrm_demo$ ll
    total 24
    drwxrwxr-x   4 craiga craiga  4096 Sep  5 16:53 ./
    drwxr-xr-x 163 craiga craiga 12288 Sep  5 16:49 ../
    drwxrwxr-x   4 craiga craiga  4096 Sep  5 16:53 recipe_kits/
    drwxrwxr-x   5 craiga craiga  4096 Sep  5 16:50 recipesystem/

The recipe execution and other commands must be added to the path. Edit your `.bashrc` or equivalent with the following exports. This is especially true since some of the commands need to be able to run each other. **Note, `/home/craiga/nrm_demo` should be changed to your installation location**.

    export PATH=**/home/craiga/nrm_demo/**recipesystem/trunk/astrodata/scripts:$PATH

The `astrodata` package needs to be in the PYTHONPATH

    export PYTHONPATH=/home/craiga/nrm_demo/recipesystem/trunk:$PYTHONPATH
    
The recipe_kits need to be added to the `RECIPEPATH` (or `ADCONFIGPATH`)

    export RECIPEPATH=/home/craiga/nrm_demo/recipe_kits
    
Source your .bashrc or otherwise execute the exports above for your own environment.
    
You can see if the system is working by using the listPrimitives.py command. 

    listPrimitives.py -i
    
This command will list transformation primitives from `kit_Novem`, which was cloned from the `recipe_kit.git` repository. Running this ensures you have the scripts directory in the PATH and `astrodata` in PYTHONPATH. If the primitives below appear, this also shows the kits are being loaded. Which should list available transformation primitives from the kit_Novem.

For a 101 Tutorial see [QUICK TUTORIAL](QUICK_TUTORIAL.md).


`listPrimitives.py` output:

```
===============================================================================
MetroBusiness 
===============================================================================
1. collapse_naics
2. convertToPercent
3. highlight
4. naics_interpret
5. sort
6. summarize_naics

        -------- Inheritance by Method Resolution Order --------

    (Pandas)
    7. columnRelate
    8. loadTables
    9. plot
    10. plotly
    11. setStorage
    12. showTables
    13. summarizeTables

    (SetRef)
    14. adaptSetType
    15. emitQAReport
    16. filterNot
    17. filterOutNot
    18. goInteractive
    19. ingest
    20. markAsIngested
    21. nativeStorage
    22. parseAsSpecial
    23. publish
    24. reduceToHeader
    25. relateData
    26. seekRelationships
    27. showContext
    28. showInputs
    29. stop
    30. writeAndDrop
    31. writeOutput
    32. writeOutputs

===============================================================================
Pandas 
===============================================================================
1. columnRelate
2. loadTables
3. plot
4. plotly
5. setStorage
6. showTables
7. summarizeTables

        -------- Inheritance by Method Resolution Order --------

    (SetRef)
    8. adaptSetType
    9. emitQAReport
    10. filterNot
    11. filterOutNot
    12. goInteractive
    13. ingest
    14. markAsIngested
    15. nativeStorage
    16. parseAsSpecial
    17. publish
    18. reduceToHeader
    19. relateData
    20. seekRelationships
    21. showContext
    22. showInputs
    23. stop
    24. writeAndDrop
    25. writeOutput
    26. writeOutputs

===============================================================================
SetRef 
===============================================================================
1. adaptSetType
2. emitQAReport
3. filterNot
4. filterOutNot
5. goInteractive
6. ingest
7. markAsIngested
8. nativeStorage
9. parseAsSpecial
10. publish
11. reduceToHeader
12. relateData
13. seekRelationships
14. showContext
15. showInputs
16. stop
17. writeAndDrop
18. writeOutput
19. writeOutputs

===============================================================================
Txt 
===============================================================================
1. parseAsDataDictionary
2. parseAsSpecial (OVERRIDES SetRef)

        -------- Inheritance by Method Resolution Order --------

    (SetRef)
    3. adaptSetType
    4. emitQAReport
    5. filterNot
    6. filterOutNot
    7. goInteractive
    8. ingest
    9. markAsIngested
    10. nativeStorage
    11. parseAsSpecial (OVERRIDDEN BY Txt)
    12. publish
    13. reduceToHeader
    14. relateData
    15. seekRelationships
    16. showContext
    17. showInputs
    18. stop
    19. writeAndDrop
    20. writeOutput
    21. writeOutputs
===============================================================================

```


<sup>*</sup> Novem is Novem, llc, of Bar Harbor, Maine at http://novem.technology/.
