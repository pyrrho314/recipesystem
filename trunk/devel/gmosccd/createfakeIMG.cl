# Scipt to take 1 noisy data frame and create faked science (with stars), flats and bias frames - MS

procedure createfakeIMG (template_file)

char    template_file    {"full_frame20_2011-05-04.fits", prompt="Test data to use to create fake 2d spectrum"}
int     num_out_files    {3,prompt="Number of output files for each type - fake spectra and flats"}
int     num_stars        {50,prompt="Number of stars per amplifier"}

begin

    int i, j, XmAx, YmAx, XmIn, YmIn, outval
    char ininfile, infile, outfile, starfile, fmt, zone, mystem, outstem
    char biasfile, flatfile, scifile
    struct udate 
    char sdate
    real bkg, mySciValue, myBiasValue, myFlatValue, noisyValue
    int masktype, num_files, l_num_stars

    #artdata
    #gemini
    #gmos
    #gemtools
    #fitsutil

    # Define the number of each type of file requested
    num_files = num_out_files
    l_num_stars = num_stars

    # Set up input name
    ininfile = template_file
    infile = mktemp("infile")
    infile = infile//".fits"

    # Define imaging keyword
    masktype = 0

    # Initiate counters
    i=1
    j=1

    # Loop of j up to the nukber of files
    for (j=1; j<=num_files; j+=1) {

        # Set up for partially random output names
        fmt = "+%S"
        zone = "-u" 
        date ( "", fmt)
        date ( "", fmt) | scan(sdate)
        outval = int(sdate) + 1 + j

        # Set multiplicative values for the faked data
        #bkg = 35300.0
        noisyValue = 49000.0
        mySciValue = (1000.0 / noisyValue ) * 1.23 * outval
        myBiasValue = (800.0 / noisyValue ) * 0.98 * outval
        myFlatValue = (4800.0 / noisyValue ) * 1.12 * outval
        #bkg = noisyValue * ( mySciValue / 1.23 ) * outval
        #bkg = mySciValue
        bkg = 1000.0

        # Copy input to a tmp file
        fxcopy ( input=ininfile, output=infile, groups="0-12" )

        # Update the header so it's imaging data
        gemhedit (infile//"[0]", "MASKTYP", masktype, "")

        # Set up a few more file nales
        mystem = "fakestars"
        udate = " "
        sdate = " "
        outstem = "ff12ampfstar"//outval
        outfile = " "
        scifile = " "
        biasfile = " "
        flatfile = " "

        # Dimensions for the locations within which the fake stars will be put
        XmAx = 512 - 80 
        YmAx = 4224 - 50 
        XmIn = 1 + 80 
        YmIn = 1 + 50 

        # Obtain time to slightly randomise the names
        fmt = "+T%H%M%S"
        zone = "-u" 
        date ( "", fmt) | scan(udate)

        # Set output names
        outfile = outstem//"-"//udate//".fits"
        scifile = "ff12ampsci"//outval//"-"//udate//".fits"
        biasfile = "ff12ampbias"//outval//"-"//udate//".fits"
        flatfile = "ff12ampflat"//outval//"-"//udate//".fits"

        # Initiate output files
        fxcopy ( input=infile, output=outfile, groups="0-12" )
        gemarith (infile,"*",myBiasValue,biasfile)
        gemarith (infile,"*",myFlatValue,flatfile)


        for ( i = 1 ; i <= 12 ; i = i + 1 ) {

            starfile = mystem//i//"_"//udate//".txt"
            print (starfile)
            # Create fake stars
            starlist (starfile, \
                l_num_stars, "", "", interactive=no, spatial="uniform", \
                xmin=XmIn., \
                xmax=XmAx., ymin=YmIn., \
                ymax=YmAx., xcenter=INDEF, ycenter=INDEF, core_radius=30., \
                base=0., \
                sseed=INDEF, luminosity="powlaw", minmag=13., maxmag=20., \
                mzero=-4., \
                power=0.6, alpha=0.74, beta=0.04, delta=0.294, mstar=1.28, \
                lseed=INDEF, nssample=100, sorder=10, nlsample=100, \
                lorder=10, \
                rbinsize=10., mbinsize=0.5, graphics="stdgraph", cursor="")

            print (infile//"["//i//"]")
            print (outfile//"["//i//"]")

            # Place fake stars into output file
            mkobjects (outfile//"["//i//"]", \
                output="", title="", ncols=512, nlines=4224, \
                header="artdata$stdheader.dat", background=bkg, 
                objects=starfile, xoffset=0.0, yoffset=0.0, \
                star="moffat", radius=6.5, beta=2.5, ar=1., pa=0., \
                distance=1., exptime=40., magzero=27., gain=1.3, rdnoise=2.3, \
                poisson=no, seed=INDEF, comments=yes)

        }

        #gemhedit (outfile//"[0]", "DETECTOR", "Hamamatsu", "Detector name")
        #gemhedit (biasfile//"[0]", "DETECTOR", "Hamamatsu", "Detector name")
        #gemhedit (flatfile//"[0]", "DETECTOR", "Hamamatsu", "Detector name")

        # Create final science frame
        gemarith (outfile,"*",mySciValue,scifile)

        # Clean up
        delete ("fakestars*", verify-, >& "dev$null")
        imdelete (infile//", "//outfile, verify-, >& "dev$null")
        #fxcopy ( input=outfile, output=scifile, groups="0-12" )
    }
END

