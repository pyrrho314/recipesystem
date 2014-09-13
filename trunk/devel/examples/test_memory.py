import numpy as np
from astrodata import AstroData
from guppy import hpy
h = hpy()

# print h.heap()
def trans_data(ad):

    for ext in ad:
        name = ext.extname()
        if name not in ['SCI','VAR','DQ']:
            continue
        print name
        data = ext.data
        print 'data accessed. (return to continue)'; 
        
        #raw_input()
        shape = (7341,7735)
        new_data = np.zeros(shape)
        print 'new data made. (return to continue)'; 
        #raw_input()
        print '\text.data id', id(ext.data)
        print '\tnew_data id', id(new_data)
        ext.data = new_data
        print '\tnew ext.data id', id(ext.data)
        del data
        print 'old data deleted. (return to continue)'; 
        #raw_input()
    ad.close()
    print 'ad closed. (return to continue)'; 
    # raw_input()


#ad0 = AstroData('../../../test_data/imcoadd/ds_mrgN20100913S0365.fits')
#ad1 = AstroData('../../../test_data/imcoadd/ds_mrgN20101002S0188.fits')
i=0
print "Press Return to start:"
raw_input()
while True:
    ad0 = AstroData('../../../test_data/imcoadd/mrgS20110228S0053.fits')
    ad1 = AstroData('../../../test_data/imcoadd/mrgS20110228S0054.fits')

    trans_data(ad0)
    trans_data(ad1)
    print "pass %d (hit return):" % i
    # raw_input()
    i += 1
    from time import sleep
    sleep(5)        
# print h.heap()

