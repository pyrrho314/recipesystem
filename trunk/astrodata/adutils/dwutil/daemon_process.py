from multiprocessing import Process, Queue
from astrodata import Lookups
from astrodata.adutils import ksutil as ks

 
class ShelfWatch(object):
    pass

ingest_sources = Lookups.compose_multi_table(
                            "*/warehouse_daemon.py",
                            "ingest_sources")

print "d_p8:", ks.dict2pretty("i_s", ingest_sources)

watch_procs = {}
if ingest_sources:
    for source in ingest_sources:
        watchproc = ShelfWatch(source)
        watch_procs[watchproc["key_path"]] = watchproc


