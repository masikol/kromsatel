# -*- encoding: utf-8 -*-

def parse_primers_lengths(primers_fpath):
    # Function parses lengths of primers stored in `primers_fpath`.
    # :param primers_fpath: path to CSV file with primers;
    # :type primers_fpath: str;
    # Returns tuple of integers where i'th integer is length
    #   of i'th primer in the input file.

    primers_len_list = list()
    min_required_cols = 2

    # We will "stretch" our primers a bit.
    # The value for `primer_len_amendment` was chosen empirically.
    primer_len_amendment = 1.25

    with open(primers_fpath, 'rt') as primers_file:
        for i, line in enumerate(primers_file):
            line_vals = tuple(map(str.strip, line.split(',')))

            if len(line_vals) < min_required_cols:
                print('Error: cannot parse line #{} in file `{}.`'.format(i+1, primers_fpath))
                print('Not enough comma-separated columns:')
                print('  {} column(s) found, {} required'.format(len(line_vals), min_required_cols))
                raise ValueError
            # end if

            primers_len_list.append(
                int(
                    len(line_vals[1]) * primer_len_amendment
                )
            )
        # end for
    # end with

    return tuple(primers_len_list)
# end def parse_primers_lengths
