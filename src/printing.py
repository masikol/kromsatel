
import sys
from time import time, strftime, gmtime

START_TIME = time() # consider time of importing as start time


def getwt():
    return strftime('%H:%M:%S', gmtime( time() - START_TIME))
# end def


def print_err(text=''):
    print(text, file=sys.stderr)
# end def
