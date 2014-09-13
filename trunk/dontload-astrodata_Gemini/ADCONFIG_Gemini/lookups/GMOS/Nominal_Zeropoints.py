nominal_zeropoints = {
    # Table of GMOS Nominal Zeropoint magnitudes
    # By filter and detector ID.
    # For new CCDs, just add them with the new detector ID
    # If we don't have per-chip values, just add the same value
    #   for each detector in the focal plane.
    #
    # The "composite" detector name is the DETID header and should here give the mean value
    # - This is used for mosaiced images that have the three detectors stiched together.
    #
    # At the moment, we don't account for the colour terms of the zeropoints
    # they are mostly < 0.1 mag,  especially longer than u, though it would be good to add
    # support for these in the future
    #
    # Nominal extinction values for CP and MK are given in a separate table 
    # at the Gemini level, and provided for by a descriptor
    #
    # Columns below are given as:
    # ("detector ID", "filter") : zeropoint
    #
    # GMOS-N original EEV detectors
    # Values from http://www.gemini.edu/sciops/instruments/gmos/calibration?q=node/10445  20111021
    ('EEV 9273-20-04', "u") : 25.47,
    ('EEV 9273-16-03', "u") : 25.47,
    ('EEV 9273-20-03', "u") : 25.47,
    ('EEV9273-16-03EEV9273-20-04EEV9273-20-03', "u") : 25.47,
    ('EEV 9273-20-04', "g") : 27.95,
    ('EEV 9273-16-03', "g") : 27.95,
    ('EEV 9273-20-03', "g") : 27.95,
    ('EEV9273-16-03EEV9273-20-04EEV9273-20-03', "g") : 27.95,
    ('EEV 9273-20-04', "r") : 28.20,
    ('EEV 9273-16-03', "r") : 28.20,
    ('EEV 9273-20-03', "r") : 28.20,
    ('EEV9273-16-03EEV9273-20-04EEV9273-20-03', "r") : 28.20,
    ('EEV 9273-20-04', "i") : 27.94,
    ('EEV 9273-16-03', "i") : 27.94,
    ('EEV 9273-20-03', "i") : 27.94,
    ('EEV9273-16-03EEV9273-20-04EEV9273-20-03', "i") : 27.94,
    ('EEV 9273-20-04', "z") : 26.78,
    ('EEV 9273-16-03', "z") : 26.78,
    ('EEV 9273-20-03', "z") : 26.78,
    ('EEV9273-16-03EEV9273-20-04EEV9273-20-03', "z") : 26.78,
    #
    # GMOS-N New (Nov 2011) deep depletion E2V CCDs.
    # Bogus values made up by PH for u-band
    ('e2v DD ML2AR CCD42-90-1-F43 10031-23-05', "u") : 10.00,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-01-03', "u") : 10.00,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-18-04', "u") : 10.00,
    ('e2v 10031-23-05,10031-01-03,10031-18-04', "u") : 10.00,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-23-05', "g") : 28.27,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-01-03', "g") : 28.27,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-18-04', "g") : 28.27,
    ('e2v 10031-23-05,10031-01-03,10031-18-04', "g") : 28.27,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-23-05', "r") : 28.37,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-01-03', "r") : 28.37,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-18-04', "r") : 28.37,
    ('e2v 10031-23-05,10031-01-03,10031-18-04', "r") : 28.37,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-23-05', "i") : 28.43,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-01-03', "i") : 28.43,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-18-04', "i") : 28.43,
    ('e2v 10031-23-05,10031-01-03,10031-18-04', "i") : 28.43,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-23-05', "z") : 27.73,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-01-03', "z") : 27.73,
    ('e2v DD ML2AR CCD42-90-1-F43 10031-18-04', "z") : 27.73,
    ('e2v 10031-23-05,10031-01-03,10031-18-04', "z") : 27.73,

    #
    # GMOS-S
    # Updated values using mean of 2013-09-02 values from website above.
    # CCD1 and CCD3 values not recoreded - changed to have the same relative
    # value to CCD2 as previous values. - MS (2013-09-04)
    ('EEV 2037-06-03', 'u'): 24.809,
    ('EEV 8194-19-04', 'u'): 24.786,
    ('EEV 8261-07-04', 'u'): 24.814,
    ('EEV2037-06-03EEV8194-19-04EEV8261-07-04', 'u'): 24.796,
    ('EEV 2037-06-03', 'g'): 28.331,
    ('EEV 8194-19-04', 'g'): 28.340,
    ('EEV 8261-07-04', 'g'): 28.371,
    ('EEV2037-06-03EEV8194-19-04EEV8261-07-04', 'g'): 28.340,
    ('EEV 2037-06-03', 'r'): 28.371,
    ('EEV 8194-19-04', 'r'): 28.361,
    ('EEV 8261-07-04', 'r'): 28.396,
    ('EEV2037-06-03EEV8194-19-04EEV8261-07-04', 'r'): 28.361,
    ('EEV 2037-06-03', 'i'): 27.941,
    ('EEV 8194-19-04', 'i'): 27.946,
    ('EEV 8261-07-04', 'i'): 27.980,
    ('EEV2037-06-03EEV8194-19-04EEV8261-07-04', 'i'): 27.946,
    ('EEV 2037-06-03', 'z'): 27.819,
    ('EEV 8194-19-04', 'z'): 26.828,
    ('EEV 8261-07-04', 'z'): 27.858,
    ('EEV2037-06-03EEV8194-19-04EEV8261-07-04', 'z'): 26.828,
}
