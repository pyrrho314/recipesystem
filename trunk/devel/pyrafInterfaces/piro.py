from astrodata import AstroData
from astrodata.RecipeManager import RecipeLibrary, ReductionContext
from astrodata import ReductionObjects # to stuff module level, should use func
from astrodata.ReductionObjects import command_clause
from astrodata.adutils import gemLog
from astrodata import Proxies

reductionObject = None
adccpid = None
prs = None

log = gemLog.getGeminiLog(logLevel=6)

def mkRO(dataset="", astrotype="", args=None, reuseRO=False):
    
    global reductionObject, adccpid
    
    if reductionObject == None or reuseRO == False:
        adccpid = Proxies.start_adcc()
        reduceServer = Proxies.ReduceServer()
        # Playing with module, should use access function
        ReductionObjects.prs = Proxies.PRSProxy.get_adcc(reduce_server=reduceServer)
        
        rl = RecipeLibrary()
        if dataset != "":
            ad = AstroData(dataset)
            ro = rl.retrieve_reduction_object(ad)
        elif astrotype != "":
            ad = None
            ro = rl.retrieve_reduction_object(astrotype = astrotype)

        # using standard command clause supplied in RecipeLibrary module
        ro.register_command_clause(command_clause)
        rc = ReductionContext()
        rc.ro = ro
        # rc.stackFile =  File("stackIndexFile", stkindfile)
        #if args:
        #    rc.add_input(ad)
        
        reductionObject = ro
    
    if reuseRO == True:
        ro = reductionObject
        rc = ro.context

    #if args:
    #    rc.addInputs(args)
    rc.update(argv)
    ro.init(rc)
    
    return reductionObject





