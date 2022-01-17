
import os
import sys
from itertools import dropwhile, takewhile


def print_help(version, last_update_data):
    print('kromsatel')
    print('Version {}. {} edition.'.format(version, last_update_data))

    readme_fpath = os.path.join(
        os.path.dirname(os.path.dirname(__file__)),
        'README.md'
    )

    with open(readme_fpath, 'rt') as readme_file:
        readle_lines = readme_file.readlines()
    # end with

    lines_after_usage = dropwhile(
        lambda line: line != '## Usage\n',
        readle_lines
    )

    lines_of_usage = takewhile(
        lambda line: line != '#### With all possible options\n',
        lines_after_usage
    )

    lines_to_rm = {'```\n',}

    lines_of_usage = filter(
        lambda line: not line in lines_to_rm,
        lines_of_usage
    )

    print()
    sys.stdout.write(''.join(lines_of_usage))
# end def print_help
