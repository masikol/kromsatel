
import os
import sys
import multiprocessing as mp

from src.printing import getwt


class Progress:

    def __init__(self, num_reads_total):
        self.NUM_READS_TOTAL = num_reads_total
        self._REPORT_DELAY = round(self.NUM_READS_TOTAL * 0.01)
        self._DEFAULT_STATUS_BAR_LEN = 40

        self.num_done_reads = mp.Value('i', 0)
        self.next_report_num = mp.Value('i', self._REPORT_DELAY)
    # end def __init__

    def get_num_done_reads(self):
        return self.num_done_reads.value
    # end def get_num_done_reads

    def get_next_report_num(self):
        return self.next_report_num.value
    # end def get_next_report_num

    def increment_done(self, increment=1):
        self.num_done_reads = self.num_done_reads + increment
    # end def increment_done

    def increment_next_report(self):
        self.next_report_num = self.next_report_num + self._REPORT_DELAY
    # end def increment_next_report

    def print_status_bar(self):
        bar_len = self._get_status_bar_len()
        curr_num_done_reads = self.get_num_done_reads()
        percent_done = round(curr_num_done_reads / self.NUM_READS_TOTAL * 100)

        sys.stdout.write(
            '\r{} - [{}] {}/{} ({}%)\n'.format(
                getwt(),
                '=' * bar_len,
                curr_num_done_reads,
                self.NUM_READS_TOTAL,
                percent_done
            )
        )
        sys.stdout.flush()
    # end def print_status_bar

    def _get_status_bar_len(self):
        try:
            bar_len = int(os.get_terminal_size().columns * 0.40)
        except OSError:
            bar_len = self._DEFAULT_STATUS_BAR_LEN
        # end try
        return bar_len
    # end _get_status_bar_len
# end class Progress
