recipesystem
============

The code in the *recipesystem* repository is released under the Mozilla Public License v. 2. The software is derived from the Gemini Observatory Python Package, sans Gemini Specific configuration packages and dependencies. Gemini released this system at ADASS 2013 in Kona, Hawaii, under a BSD license and subsequently, and it is expected to arrive on GitHub or some other public git repository.

This version is released under MPL 2. 

At the present time, Novem LLC is generalizing the system for commercial use and using the system to automate its clients' particular data processing pipelines using the recipe model. Key development goals for the immediate future:

* relocation of AstroData source code into the astrodata_FITS configuration package
* general purpose documentation (available documentation is Astronomy Oriented)
* automatic cloud-based scaling, job cooperation
* general purpose HTML5/Javascript interfaces for system control and monitoring
* generic library of numpy and pandas primitives
* expansion of the general purpose Set Reference Data Type
* creation of a project roadmap

Install
-------

I'm creating a directory `nrm_demo` and executing the following there.

First get the instrastructure in the Novem Recipe System from git.

    git clone https://github.com/pyrrho314/recipesystem.git

We'll add a script directory to the path, and a python directory to PYTHONPATH below, but first, get the standard recipe kit, kit_Novem.

    git clone https://github.com/pyrrho314/recipe_kits.git
    
This gives us:
    craiga@pleiades:~/nrm_demo$ ll
    total 24
    drwxrwxr-x   4 craiga craiga  4096 Sep  5 16:53 ./
    drwxr-xr-x 163 craiga craiga 12288 Sep  5 16:49 ../
    drwxrwxr-x   4 craiga craiga  4096 Sep  5 16:53 recipe_kits/
    drwxrwxr-x   5 craiga craiga  4096 Sep  5 16:50 recipesystem/

The recipe execution and other commands must be added to the path. This is especially true since some of the commands need to be able to run each other. Note, `/home/craiga/nrm_demo` should be changed to your installation location.

    export PATH=/home/craiga/nrm_demo/recipesystem/trunk/astrodata/scripts:$PATH

The `astrodata` package needs to be in the PYTHONPATH

    export PYTHONPATH=/home/craiga/nrm_demo/recipesystem/trunk:$PYTHONPATH
    
The recipe_kits need to be added to the `RECIPEPATH` (or `ADCONFIGPATH`)

    export RECIPEPATH=/home/craiga/nrm_demo/recipe_kits
    
You can see if the system is working by typeing:

    listPrimitives.py
    
This command will list transformation primitives from `kit_Novem`, which was cloned from the `recipe_kit.git` repository. Running this ensures you have the scripts directory in the PATH and `astrodata` in PYTHONPATH. If the primitives below appear, this also shows the kits are being loaded. Which should list available transformation primitives from the kit_Novem:

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
