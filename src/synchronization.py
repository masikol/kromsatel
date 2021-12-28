
import multiprocessing as mp


num_done_reads  = mp.Value('i', 0)
next_report_num = mp.Value('i', 0)

output_lock        = mp.Lock()
print_lock         = mp.Lock()
status_update_lock = mp.Lock()
