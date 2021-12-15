
import multiprocessing as mp

class Synchronizer:

    def __init__(self):
        self.output_lock = mp.Lock()
        self.print_lock = mp.Lock()
        self.status_update_lock = mp.Lock()
    # end def __init__
# end class Synchronizer
