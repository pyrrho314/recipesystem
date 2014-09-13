import sys, os
import logging
import datetime

def stringify(args):
    return " ".join(map(str, args))

def debug(*args):
    log = logging.getLogger("GACQ")
    log.debug(stringify(args))

def info(*args):
    log = logging.getLogger("GACQ")
    log.info(stringify(args))

def error(*args):
    log = logging.getLogger("GACQ")
    log.error(stringify(args))

def warning(*args):
    log = logging.getLogger("GACQ")
    log.warning(stringify(args))

def is_debug_mode():
    log = logging.getLogger("GACQ")
    return log.getEffectiveLevel() == logging.DEBUG

def setup_logging(gacqlog, verbose):
    loglevel = logging.INFO
    if verbose:
        loglevel = logging.DEBUG

    log = logging.getLogger("GACQ")
    log.setLevel(loglevel)

    if gacqlog is None or gacqlog == "":
        base, ext = os.path.splitext(os.path.basename(sys.argv[0]))
        datestring = datetime.datetime.utcnow().strftime("%Y-%m-%d")
        gacqlog = base + ".log." + datestring # Log file

    # log file 
    logfile = logging.FileHandler(gacqlog)
    logfile.setLevel(loglevel)
    formatter = logging.Formatter(fmt="%(asctime)s %(name)-12s %(message)s",
                                  datefmt="%Y-%m-%d %H:%M:%S")
    logfile.setFormatter(formatter)
    log.addHandler(logfile)

    # console
    console = logging.StreamHandler()
    console.setLevel(loglevel)
    formatter = logging.Formatter(fmt="%(name)-12s %(message)s")
    console.setFormatter(formatter)
    log.addHandler(console)

    debug("...log = ", gacqlog)

    return logfile, console
