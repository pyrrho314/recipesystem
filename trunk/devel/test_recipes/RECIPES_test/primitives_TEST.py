from primitives_GEMINI import GEMINIPrimitives
\
from astrodata.ReductionObjects import PrimitiveSet

class TestPrimitives (GEMINIPrimitives):
    def stepOne(self, rc):
        print "stepOne"
        print rc["stepOne"]
        
        print rc["arg"]
        yield rc
        
    def stepTwo(self, rc):
        print "stepTwo"
        print rc["stepTwo"]
        print rc["arg"]
        yield rc
    def stepThree(self, rc):
        print "stepThree"
        print rc["stepThree"]
        print rc["arg"]
        yield rc  
        
    def stepCaller(self, rc):
        print "stepCalling running stepTwo, then called"
        rc.run("stepTwo(arg=8) #love ya")
        rc.run("called(arg=888)")    
        if True:
            rc.run("""
            stepOne(arg=1)
            stepTwo(arg=2)
            showParams(arg=True)
                """)
        yield rc
