.. FLUX_CAL  Documentation master file

=====================
FLUXCAL Documentation
=====================

   Fluxcal consists of a set of scripts that lead to determine the Zero Point correction
   value from image fields in a FITS file. The steps require to achieve this are:

   - Detects object sources. 
       Creates a FITS table 'OBJCAT' with x,y,ra,dec,flux for each source in the field.
   - Adding Reference Catalog sources.
       Creates a FITS table 'REFCAT' with ra,dec,x,y,mag for each standard object
       that falls in the image field.
   - Correlate Objects and Reference sources.
       Finds reference stars that matches object positions in the field within
       a given radious.
   - Calculate Zero point.
       Computes photometry for those objects found in the correlation.

.. toctree::
   :maxdepth: 1

   running_fc

   scripts

.. **Detect Source Algoritm**
    --------------------------------

    .. automodule:: detSources
    .. autoclass:: detSources
        :members:


Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`

