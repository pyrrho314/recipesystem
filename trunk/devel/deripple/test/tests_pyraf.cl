# 1) start pyraf
# 2) # Set deripplepath to where the deripple.py module is
#    # where the module is located
#    set deripplepath=<some_path>/deripple/
# 3) # Define the pyraf task
#    task deripple = deripplepath$deripple.cl
#
# 4) # Load the task into the pyraf environment
#    deripple
#
# 5) # Run the task

deripple(infile='gn20120419_104_sci.fits',outfile='out_104.fits',
               start=550 ,end=610 ,period=2)

deripple(infile='gn20120419_84_sci.fits',outfile='out_84.fits',
               start=550 ,end=610 ,period=2)

