
.. _Detect-Sources:

Fluxcal Test script
-------------------

.. automodule:: fluxcal
    :members:

Detect Sources
---------------

  This is the first script that you need to run to be able to get 
  result in CalculateZeropoint. It uses the 'astrodata' package to 
  handle input/output and a detection algorithm that follows daofind.
  Please see :ref:`Det-Sources`.

.. automodule:: detectSources
    :members:

Add Reference Catalog
-----------------------
.. automodule:: addReferenceCatalogs
.. autoclass:: AddReferenceCatalogs
    :members:

Correlate Object and reference positions
-----------------------------------------
.. automodule:: correlateWithReferenceCatalogs
.. autoclass:: CorrelateWithReferenceCatalogs
    :members:

Calculate Zeropoint, Image quality
-----------------------------------------

.. automodule:: calculateZeropoint
    :members:

.. _Det-Sources:

DetSources
----------
.. automodule:: detSources
    :members: detSources

Primitives
----------
.. automodule:: primitives_fluxcal_GMOS_IMAGE
    :members: GMOS_IMAGE_fluxcal_Primitives
