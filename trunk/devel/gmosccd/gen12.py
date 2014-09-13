#!/usr/bin/env python

import os, sys
from optparse import OptionParser
from time import sleep
from pyraf import iraf
from pyraf.iraf import stsdas, tables
from pyraf.iraf import gemini
from pyraf.iraf import gmos
from pyraf.iraf import images
from pyraf.iraf import gemtools
import pyfits as pf
import datetime
# Parse command line
parser = OptionParser()
parser.set_description("""Generate 12ext fits""")

parser.add_option( '-d', '--display', action='store_true', dest='display',
                   default=False, help='Displays end image using gdisplay' )
parser.add_option( '-m', '--metadata', dest='metadata',
                   default=False, help='Will display phu according to given parameter, 0=phu 1=1st sci ext' )
( options, args ) = parser.parse_args()

# set up keyword arrays ( without binning dependency)
detsec=[   '1:512',  '513:1024','1025:1536','1537:2048',
        '2049:2560','2561:3072','3073:3584','3585:4096',
        '4097:4608','4609:5120','5121:5632','5633:6144']
ccdName=['EEVccd_001', 'EEVccd_002', 'EEVccd_003']
gain=[2.03999996185303, 2.3199999332428, 2.19000005722046]
rdNoise=[3.5, 3.29999995231628, 3]

# begin loop through images
for arg in args:
    a3 = pf.open( arg )

    # sort out binning
    oneByOne=False
    twoByTwo=False
    twoByOne=False
    if a3[1].header['CCDSUM'] == '1 1':
        oneByOne = True
    elif a3[1].header['CCDSUM'] == '2 2':
        twoByTwo = True
    elif a3[1].header['CCDSUM'] == '2 1':
        twoByOne = True
    else:
        print 'Only will accept binning (CCDSUM) 1 1 or 2 2 or 2 1'
        sys.exit()
    os.system('rm slice*')
    trimmed=False
    try:
        if a3[0].header['TRIMMED'] == 'yes':
            print "Image was trimmed..."
            trimmed=True
    except:
        print "Image NOT trimmed..."
        pass

    # set up slicing parameters, copy overscan regions 
    # and extension keyword arrays (with binning dependency)
    #os.system('rm slice*')
    if trimmed:
        if oneByOne:
            print 'Image has been trimmed, binning = 1x1.....'
            xcopy1_2 =[ '33:544',  '545:1056','1057:1568','1569:2080' ]
            xcopy_3  =[  '1:512',  '513:1024','1025:1536','1537:2048' ]
            datasec=['[33:544,1:4608]','[1:512,1:4608]']
            biassec=['[1:32,1:4608]','[513:544,1:4608]']
        elif twoByTwo:
            print 'Image has been trimmed, binning = 2x2.....'
            #xcopy1_2=[   '33:288', '289:544', '545:800', '801:1056']
            xcopy_3 =[    '1:256', '257:512', '513:768', '769:1024']
            xcopy_1_2 = xcopy_3
            datasec=['[1:256,1:2304]','[1:256,1:2304]']
            biassec=['[1:32,1:2304]','[257:288,1:2304]']
        else:
            print 'Binning = 2x1.....'
            print 'Image has been trimmed, binning = 2x1.....'
            #xcopy1_2=[   '33:288', '289:544', '545:800', '801:1056']
            xcopy_3 =[    '1:256', '257:512', '513:768', '769:1024']
            xcopy_1_2 = xcopy_3
            datasec=['[1:256,1:4608]','[1:256,1:4608]']
            biassec=['[1:32,1:4608]','[257:288,1:4608]']
    else:
        images.imcopy( arg+'[1][1:32,*]'   ,'sliceOverscan1' )
        images.imcopy( arg+'[2][1:32,*]'   ,'sliceOverscan2' )
        if oneByOne:
            print 'Binning = 1x1.....'
            xcopy1_2 =[ '33:544',  '545:1056','1057:1568','1569:2080' ]
            xcopy_3  =[  '1:512',  '513:1024','1025:1536','1537:2048' ]
            datasec=['[33:544,1:4608]','[1:512,1:4608]']
            biassec=['[1:32,1:4608]','[513:544,1:4608]']
            images.imcopy( arg+'[3][2049:2080,*]','sliceOverscan3' )
                
        elif twoByTwo:
            print 'Binning = 2x2.....'
            xcopy1_2=[   '33:288', '289:544', '545:800', '801:1056']
            xcopy_3 =[    '1:256', '257:512', '513:768', '769:1024']
            datasec=['[33:288,1:2304]','[1:256,1:2304]']
            biassec=['[1:32,1:2304]','[257:288,1:2304]']
            images.imcopy( arg+'[3][1025:1056,*]','sliceOverscan3' )
        else:
            print 'Binning = 2x1.....'
            xcopy1_2=[   '33:288', '289:544', '545:800', '801:1056']
            xcopy_3 =[    '1:256', '257:512', '513:768', '769:1024']
            datasec=['[33:288,1:4608]','[1:256,1:4608]']
            biassec=['[1:32,1:4608]','[257:288,1:4608]']
            images.imcopy( arg+'[3][1025:1056,*]','sliceOverscan3' )
            
    # slice raw image using 'imcopy' then join them with overscan using 'imjoin'
    print '\nCreating Ext12_'+arg,'.....'
    count=0
    for i in range(1,4):
        ext='['+str(i)+']['
        over='sliceOverscan'+str(i)+'.fits'
        for j in range(4):
            count+=1
            pre_sl='slicePre'+str(count)
            sl='slice'+str(count)
            if trimmed:            #no need to join if image alread trimmed
                    images.imcopy( arg+ext+xcopy_3[j]+',*]',sl )
            else:
                if i==3:           #third chip requires different cut because overscan on right
                    images.imcopy( arg+ext+xcopy_3[j]+',*]',pre_sl )
                else:
                    images.imcopy( arg+ext+xcopy1_2[j]+',*]',pre_sl )
                if not count%2:    #using overscan right then left convention
                    images.imjoin( input=over+','+pre_sl+'.fits',output=sl, join_dimension=1 )
                else:
                    images.imjoin( input=pre_sl+'.fits,'+over,output=sl, join_dimension=1 )
    
    # create the multi-extension fits file by running 'wmef' on the slices (w/overscan)
    print 'wmef....'
    newarg = 'Ext12_'+arg
    os.system('rm '+newarg)
    gemtools.wmef(input='slice1.fits,slice2.fits,slice3.fits,slice4.fits,slice5.fits,slice6.fits,slice7.fits,slice8.fits,slice9.fits,slice10.fits,slice11.fits,slice12.fits',output=newarg, phu=arg)
    os.system('rm slice*')

    # update header keywords in 12-ext image
    print 'updating headers...'
    a12=pf.open( newarg,mode='update' )
    a12[0].header.update( 'NAMPS'   , 4,                           'Number of amplifiers')
    a12[0].header.update( 'RELEASE' , '2002-05-07','End of proprietary period YYYY-MM-DD')
    a12[0].header.update( 'DETECTOR',                    'Hamamatsu', 'Amplifier name(s)')
    a12[0].header.update( 'DETID',  'EEVccd_001EEVccd_002EEVccd_003', 'All chip names')
    a12[0].header.update( 'DATE',datetime.datetime.utcnow().isoformat(), 'Date FITS file was generated')
    if trimmed:
        a12[0].header.update( 'NSCIEXT', 12 , 'Number of science extensions')
        
    count = 0
    for i in range(12):
        ampname=ccdName[count]+', amp'+str(i+1)
        a12[i+1].header.update( 'CCDNAME',            ccdName[count], 'CCD chip name(s)' )
        a12[i+1].header.update( 'AMPNAME',                   ampname, 'Amplifier name(s)')
        if trimmed:
            a12[i+1].header.update( 'EXTNAME',                     'SCI',    'Extension name')
            a12[i+1].header.update( 'EXTVER',                        i+1, 'Extension version')
        
        if not i%2:
            a12[i+1].header.update( 'DATASEC',                          datasec[1], 'CCD section(s)'   )
            a12[i+1].header.update( 'BIASSEC',                          biassec[1], 'CCD section(s)'   )
            if trimmed:
                if twoByTwo:
                    a12[i+1].header.update( 'TRIMSEC', '[1:256,1:2304]', ''   )
                if twoByOne:
                    a12[i+1].header.update( 'TRIMSEC', '[1:256,1:4608]', ''   )
        else:
            if trimmed:
                if twoByTwo:
                    a12[i+1].header.update( 'TRIMSEC', '[33:288,1:2304]', ''   )
                if twoByOne:
                    a12[i+1].header.update( 'TRIMSEC', '[33:288,1:4608]', ''   )
            a12[i+1].header.update( 'DATASEC',                          datasec[0], 'CCD section(s)'   )
            a12[i+1].header.update( 'BIASSEC',                          biassec[0], 'CCD section(s)'   )
        if not trimmed:
            a12[i+1].header.update( 'GAIN   ',                     gain[count] , 'Amplifier gain'   )
            a12[i+1].header.update( 'RDNOISE',                   rdNoise[count], 'Readout noise'    )
        a12[i+1].header.update( 'CCDSEC' ,         '['+detsec[i%4]+',1:4608]', 'CCD section(s)'   )
        a12[i+1].header.update( 'DETSEC' ,           '['+detsec[i]+',1:4608]', 'CCD section(s)'   )
        if i%4==3:
            count+=1
        
    a12.flush()
    a12.info()
    if options.metadata:
        print a12[int( options.metadata )].header
    a12.close()
    a3.close()
    if options.display:
        gmos.gdisplay( newarg,frame=1,fl_paste='yes' )
        
    print 'done'

#___________________________________________________________________________eof

