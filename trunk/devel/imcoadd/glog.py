import time

# These are the loglevel arguments and the corresponding strings that
# will actually be written to the log file.
logleveldict = {"engineering": "ENG", "science": "SCI", "status": "STAT",
                "task": "TSK", "visual": "VIS", "none": "NONE"}

# These are the GEMLOG error numbers and corresponding error messages.
gemerrmsg = { 99: "Internal error",
             100: "Error opening file",
             101: "Unable to access file",
             102: "File already exists",
             120: "Using default value",
             121: "Input error",
             122: "Unrecognized option",
             123: "Wrong image format",
             131: "Keyword not found",
             132: "Error in header content"
            }

class GLog:

    def __init__ (self, logfile, curtask=None, curpack=None, paramstr=None,
                  fl_append=True, glogpars=None, verbose=True):

        # curpack and glogpars are currently not used

        if curtask is not None:
            self._curtask = curtask.upper()
        else:
            self._curtask = ""

        self._fl_success = True

        # self.verbose is the one used to control printing to standard output;
        # this can be overridden if verbose is specified in glogprint or
        # glogclose.
        # self._verbose is taken from the argument to __init__, and this will
        # be copied to self.verbose in glogprint and glogclose if the
        # verbose argument to those methods was not specified explicitly.
        self._verbose = verbose
        self.verbose = verbose

        if fl_append:
            mode = "a"
        else:
            mode = "w"
        self._fd = None
        self._fd = open (logfile, mode=mode)

        self.printBegin()

        if paramstr is not None:
            self.writeParam (paramstr)

        self._fd.flush()

    def glogprint (self, str="", loglevel="status",
                   type="string", fork="forward", child="",
                   vistype="empty", errno=0, verbose=None):
        """Write a message to the log file."""

        if verbose is None:
            self.verbose = self._verbose   # use the value specified initially
        else:
            self.verbose = verbose         # override the initial value

        logstr = logleveldict[loglevel]
        if type == "string":
            self.writeStr (logstr, str)
        elif type == "file" or type == "list":
            self.writeFileContents (logstr, type, str)
        elif type == "error" or type == "warning":
            self.writeErrorMsg (logstr, type, errno, str)
        elif type == "fork":
            self.writeForkMsg (logstr, type, fork, child)
        elif type == "visual":
            self.writeVisual (logstr, vistype)
        else:
            raise RuntimeError, \
                "error in a call to glogprint:  type = %s" % type

        self._fd.flush()

    def glogclose (self, fl_success=None, verbose=None):
        """Write the final messages and close the log file."""

        if self._fd is None:            # already closed?
            return

        if fl_success is not None:
            self._fl_success = fl_success

        if verbose is None:
            self.verbose = self._verbose
        else:
            self.verbose = verbose

        self.printEnd()

        self._fd.close()
        self._fd = None

    def printBegin (self):
        """Write the initial (BOE) messages."""

        (date_time, time_date) = self.getDateTime()

        self.writeStr ("BOE", date_time)

        message = "Log opened at [%s]" % time_date
        self.writeStr ("STAT", message)

        self.writeStr ("VIS", "")

    def printEnd (self):
        """Write the final (EOE) messages."""

        (date_time, time_date) = self.getDateTime()

        if self._fl_success:
            status = "SUCCESS"
        else:
            status = "FAILURE"
        message = "Exit status: %s" % status
        self.writeStr ("STAT", message)

        message = "Log closed at [%s]" % time_date
        self.writeStr ("STAT", message)

        self.writeStr ("EOE", date_time)
        self.writeStr ("VIS")

    def writeFileContents (self, loglevel, type, str):
        """Write the contents of a text file to the log file."""

        if type == "file":
            fd = open (str)             # str is the file name
            lines = fd.readlines()
            fd.close()
        else:
            lines = str                 # str is a list of strings

        for line in lines:
            line = line.rstrip()
            self.writeStr (loglevel, line)

    def writeErrorMsg (self, loglevel, type, errno, str):
        """Write an error message to the log file."""

        self.verbose = True
        type = type.upper()

        (date_time, time_date) = self.getDateTime()
        message = "%s: %d at [%s]" % (type, errno, time_date)
        self.writeStr (loglevel, message)

        if errno != 0:
            errmess = self.getErrorMsg (errno)
            if errmess:
                message = "%s: %d %s" % (type, errno, errmess)
                self.writeStr (loglevel, message)
        if str:
            message = "%s: %d %s" % (type, errno, str)
            self.writeStr (loglevel, message)

    def getErrorMsg (self, errno):
        """Find the error message corresponding to an error number."""

        if gemerrmsg.has_key (errno):
            return gemerrmsg[errno]
        else:
            return ""

    def writeForkMsg (self, loglevel, type, fork, child):
        """Write a message regarding a function call."""

        if fork == "forward":
            message = "-- Forking to %s ..." % child.upper()
        elif fork == "backward":
            message = "-- Returning to %s ..." % self._curtask
        else:
            raise RuntimeError, \
                "error in a call to glogprint:  fork = %s" % fork
        self.writeStr (loglevel, message)

    def writeVisual (self, loglevel, vistype):
        """Write text to serve as a visual marker."""

        if vistype == "empty":
            message = ""
        elif vistype == "shortdash":
            message = "--------------------"
        elif vistype == "longdash":
            message = \
            "-----------------------------------------------------------"
        else:
            raise RuntimeError, \
                "error in a call to glogprint:  vistype = %s" % vistype
        self.writeStr (loglevel, message)

    def getDateTime (self):
        """Get the date and time strings (in two formats)."""

        local_time = time.localtime (time.time())
        date_time = time.strftime ("%Y-%m-%dT%H:%M:%S", local_time)
        time_date = time.strftime ("%a %H:%M:%S %d-%b-%Y", local_time)
        return (date_time, time_date)

    def writeStr (self, loglevel="", message=""):
        """Write a string to the log file and optionally to standard output."""

        str = "%-4s %s" % (loglevel, self._curtask)
        if message:
            str = str + " " + message

        self._fd.write (str + "\n")
        if self.verbose:
            if len (str) < 6:
                print ""
            else:
                print str[5:]

    def writeParam (self, paramstr):
        """Write a set of parameters to the log file."""

        self.glogprint (loglevel="visual", type="visual", vistype="longdash")
        self.writeStr ("TSK", "Input Parameters:")
        for line in paramstr:
            words = line.split ("=")
            message = "     %-11s=" % words[0].strip()
            if len (words) > 1:
                value = words[1].strip()
                if value:
                    message = message + " %s" % value
            self.writeStr ("TSK", message)
        self.glogprint (loglevel="visual", type="visual", vistype="longdash")
