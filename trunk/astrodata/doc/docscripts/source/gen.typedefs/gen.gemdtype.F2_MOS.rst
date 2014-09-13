
F2_MOS Classification Source
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

.. toctree::
    :numbered:
    :maxdepth: 0
     
Classification
    F2_MOS

Source Location 
    ADCONFIG_Gemini/classifications/types/F2/gemdtype.F2_MOS.py

.. code-block:: python
    :linenos:

    class F2_MOS(DataClassification):
        name="F2_MOS"
        usage = """
            Applies to all MOS datasets from the FLAMINGOS-2 instrument
            """
        parent = "F2_SPECT"
        requirement = AND ([  ISCLASS("F2_SPECT"),
                              PHU(OBSTYPE="OBJECT"),
                              OR([  PHU(DECKER="mos"),
                                    PHU(DCKERPOS="mos"),
                                    PHU(MOSPOS="mos.?")  ])  ])
    
    newtypes.append(F2_MOS())



