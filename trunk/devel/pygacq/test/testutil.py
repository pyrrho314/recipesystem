import os, sys
from cStringIO import StringIO

def get_data_file_name(fname):
    directories_to_search = []

    env_var_name = "GACQ_TEST_DATA_DIR"
    if env_var_name in os.environ:
        directories_to_search.append(os.environ[env_var_name])

    directories_to_search.append(os.path.join(os.path.dirname(__file__), "data"))
    directories_to_search.append("/rtfperm/rtftests/gemaux_python/gempylocal/orig")

    for dname in directories_to_search:
        fullpath = os.path.join(dname, fname)
        if os.path.exists(fullpath):
            return fullpath
        
    raise ValueError("Unable to locate '%s' in the following directories: %r. Either place it in one of those directories or specify where to find it with the '%s' environment variable" %
                     (fname, directories_to_search, env_var_name))

def assert_almost_equals(val1, val2, places=7, msg=None):
    import nose.tools
    try:
        iter(val1)
        try:
            iter(val2)
        except TypeError:
            assert False, "%s is iterable but %s is not" % (repr(val1), repr(val2))

        assert len(val1) == len(val2)
        if isinstance(val1, dict):
            if not isinstance(val2, dict):
                assert False, "%s is a dict but %s is not" % (repr(val1), repr(val2))
                
            for item1, item2 in zip(val1.items(), val2.items()):
                assert_almost_equals(item1, item2, places, msg)

        for v1, v2 in zip(val1, val2):
            assert_almost_equals(v1, v2, places, msg)
                
    except TypeError:
        nose.tools.assert_almost_equals(val1, val2, places, msg)

def assert_tolerance(val1, val2, tolerance=0.001, msg=None):
    try:
        iter(val1)
        try:
            iter(val2)
        except TypeError:
            assert False, "%s is iterable but %s is not" % (repr(val1), repr(val2))

        assert len(val1) == len(val2)
        if isinstance(val1, dict):
            if not isinstance(val2, dict):
                assert False, "%s is a dict but %s is not" % (repr(val1), repr(val2))
                
            for item1, item2 in zip(val1.items(), val2.items()):
                assert_tolerance(item1, item2, tolerance, msg)

        for v1, v2 in zip(val1, val2):
            assert_tolerance(v1, v2, tolerance, msg)
                
    except TypeError:
        diff = abs(val1 - val2)
        assert diff < tolerance, "|%f - %f| = %f > %f" % (val1, val2, diff, tolerance)

def setup_io(stdin=""):
    if stdin:
        stdin += "\n"
    stdin += "no"
    def setup_func():
        sys.stdin = StringIO(stdin) 
        sys.stderr = StringIO()
    return setup_func

def get_offsets():
    stdout = sys.stderr.getvalue()
    lines = stdout.split("\n")

    for line in lines:
        marker = "OFFSETS (arcsec):"
        if marker in line:
            parts = line.split(marker)
            assert len(parts) == 2

            parts = parts[1].split()

            rotmarker = "ROTATION"
            if rotmarker in line:
                assert len(parts) == 10
                return float(parts[2]), float(parts[5]), float(parts[9])
            else:
                assert len(parts) == 6
                return float(parts[2]), float(parts[5])

    assert False, "Could not find offsets in GACQ output to stderr:\n" + stdout
    
