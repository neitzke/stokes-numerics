from __future__ import absolute_import
import logging

def logconfig(console=True,filename='/tmp/harmonic.log'):
    formatter = logging.Formatter('[%(asctime)s] [%(processName)s] [%(name)s] [%(levelname)s] %(message)s', datefmt='%m-%d %H:%M')

    root = logging.getLogger('')
    root.setLevel(logging.INFO)
    
    if filename:
        # define a Handler which writes some messages to a log file
        logfilehandler = logging.FileHandler(filename)
        logfilehandler.setFormatter(formatter)

        # add the handler to the root logger
        logging.getLogger('').addHandler(logfilehandler)

    if console:
        # define a Handler which writes some messages to sys.stderr
        consolehandler = logging.StreamHandler()
        consolehandler.setFormatter(formatter)

        # add this handler to the root logger
        logging.getLogger('').addHandler(consolehandler)

# to change log level globally, use eg logconfig.loglevel(logging.WARN)
# to change level for an individual module, eg logconfig.loglevel(logging.DEBUG, "framedata")
def loglevel(level, name = ""):
    logging.getLogger(name).setLevel(level)
