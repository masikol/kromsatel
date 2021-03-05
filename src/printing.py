# -*- coding: utf-8 -*-

from time import time, strftime, gmtime

START_TIME = time() # consider time of importing as start time


def getwt():
    # Function (get work time) returns time HH:MM:SS that has passed from start_time.
    return strftime('%H:%M:%S', gmtime( time() - START_TIME))
# end def getwt
