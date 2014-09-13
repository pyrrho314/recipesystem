import os

def local_fsc_patch(fsc):
    fsc.using_sqlite = True
    fsc.using_apache = False
    fsc.fsc_localmode = True
    fsc.using_cadc = False
    fsc.storage_root = dpath = os.path.abspath(".")
    fsc.das_calproc_path = dpath 
    fsc.fits_servername = "local"
    fsc.fits_system_status = "development"
    fsc.email_errors_to = "callen@gemini.edu"
    fsc.fits_dbname = "fitsdata.sql"
    fsc.fits_database = "sqlite:///%s/%s" % (dpath, fsc.fits_dbname) 
    fsc.fits_db_backup_dir = "%s/backups" % dpath
    fsc.fits_lockfile_dir = "%s/locks" % dpath
    fsc.fits_log_dir = "%s/logs/" % dpath
    fsc.fits_tape_scratchdir = "%s/tapescratch" % dpath

    from astrodata import Descriptors
    os.environ.update({"FITSSTORAGECONFIG_LOCALMODE":"True"})
