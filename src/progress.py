
import os
import sys

from src.printing import getwt
import src.synchronization as synchron


class Progress:

    def __init__(self, num_reads_total):
        self.NUM_READS_TOTAL = num_reads_total
        self._REPORT_DELAY = max(
            1,
            round(self.NUM_READS_TOTAL * 0.01)
        )
        self._DEFAULT_STATUS_BAR_LEN = 40

        synchron.next_report_num.value = self._REPORT_DELAY
    # end def


    def get_num_done_reads(self):
        return synchron.num_done_reads.value
    # end def


    def get_next_report_num(self):
        return synchron.next_report_num.value
    # end def


    def increment_done(self, increment=1):
        synchron.num_done_reads.value = synchron.num_done_reads.value + increment
    # end def


    def increment_next_report(self):
        while synchron.num_done_reads.value >= synchron.next_report_num.value:
            synchron.next_report_num.value = synchron.next_report_num.value + self._REPORT_DELAY
        # end while
    # end def


    def print_status_bar(self):

        curr_num_done_reads = self.get_num_done_reads()

        bar_len = self._get_status_bar_len()
        ratio_done = curr_num_done_reads / self.NUM_READS_TOTAL
        percent_done = ratio_done * 100
        progress_line_len = round(bar_len * ratio_done)

        print_arrow = progress_line_len != bar_len
        if print_arrow:
            arrow = '>'
        else:
            arrow = ''
            progress_line_len += 1
        # end if

        sys.stdout.write(
            '\r{} - [{}{}{}] {}/{} ({}%)'.format(
                getwt(),
                '=' * progress_line_len,
                arrow,
                ' ' * (bar_len - progress_line_len),
                curr_num_done_reads,
                self.NUM_READS_TOTAL,
                round(percent_done)
            )
        )
        sys.stdout.flush()
    # end def


    def _get_status_bar_len(self):
        try:
            bar_len = int(os.get_terminal_size().columns * 0.40)
        except OSError:
            bar_len = self._DEFAULT_STATUS_BAR_LEN
        # end try
        return bar_len
    # end def
# end class
