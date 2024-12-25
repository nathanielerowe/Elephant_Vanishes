import profit
import inspect
import os

# Loggers

# Log stuff the way PROfit tends to in c++ world
def pyPROlog(level, s):

    # Callers should be two levels up the stack, 
    # since this should be called by one of the PROlog* functions 
    caller = inspect.stack()[2]
    # The name is in the fourth element
    caller_name = caller[3]

    # Also get the script. Make sure we have just the file name.
    script_name = os.path.basename(inspect.stack()[-1][1])

    profit.PROlog(level, "%s::%s || %s" % (script_name, caller_name, s))

def PROlogINFO(s):
    return pyPROlog(profit.globals.LOG_INFO, s)

def PROlogDEBUG(s):
    return pyPROlog(profit.globals.LOG_DEBUG, s)

def PROlogWARNING(s):
    return pyPROlog(profit.globals.LOG_WARNING, s)

def PROlogERROR(s):
    return pyPROlog(profit.globals.LOG_ERROR, s)

def PROlogCRITICAL(s):
    return pyPROlog(profit.globals.LOG_CRITICAL, s)
