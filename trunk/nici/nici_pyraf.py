from pyraf import iraf
import nici
from nici.niciTools import getFileList
from nici import ncmkflats,ncprepare,ncscience,ncqlook
import imp

def _ncmkflats (inputs, idir, odir, sigma, clobber, suffix, logfile, verbose):

    list = getFileList(inputs)
    clobber = (clobber == 'yes')
    verbose = (verbose == 'yes')

    ncmkflats.ncmkflats (list, idir, odir, sigma, clobber, suffix, logfile, verbose)
 
parfile = iraf.osfn('nicipath$ncmkflats.par')
t = iraf.IrafTaskFactory(taskname='ncmkflats', value=parfile, function=_ncmkflats)

def _ncprepare(inputs, oprefix, idir, odir, fdir, fsuffix, dobadpix, clobber, logfile, verbose):

    list = getFileList(inputs)
    dobadpix = (dobadpix == 'yes')
    clobber = (clobber == 'yes')
    verbose = (verbose == 'yes')

    ncprepare.ncprepare(list, oprefix, idir, odir, fdir, fsuffix, dobadpix, clobber, logfile, verbose)
 
parfile = iraf.osfn('nicipath$ncprepare.par')
t = iraf.IrafTaskFactory(taskname='ncprepare', value=parfile, function=_ncprepare)

def _ncscience (inputs, idir, odir, central, suffix, bsize, mdfw, clobber, logfile, verbose):

    list = getFileList(inputs)
    clobber = (clobber == 'yes')
    central = (central == 'yes')
    verbose = (verbose == 'yes')

    ncscience.ncscience (list, idir, odir, central, suffix, bsize, mdfw, clobber, logfile, verbose)
 
parfile = iraf.osfn('nicipath$ncscience.par')
t = iraf.IrafTaskFactory(taskname='ncscience', value=parfile, function=_ncscience)

def _ncqlook (inputs, idir, odir, log,lists,saturate,nodisplay,full,port):

    saturate = 5000
    log = True
    lists = True
    nodisplay = (nodisplay == 'yes')
    full = (full == 'yes')
    port = 5137
    ncqlook.ncqlook (inputs, idir, odir, log,lists,saturate, nodisplay, full, port)
 
parfile = iraf.osfn('nicipath$ncqlook.par')
t = iraf.IrafTaskFactory(taskname='ncqlook', value=parfile, function=_ncqlook)
