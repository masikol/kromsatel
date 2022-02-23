#!/usr/bin/env python3

__version__ = '2.0.b'
#                       YYYY-mm-dd
__last_update_date__ = '2022-02-23'
# __author__ = 'Maxim Sikolenko'


# Check python interpreter version

import sys

if sys.version_info.major < 3:
    print(
        '\nYour python interpreter version is ' + '%d.%d' \
            % (sys.version_info.major, sys.version_info.minor)
    )
    print('   Please, use Python 3.\a')
    # In python 2 'raw_input' does the same thing as 'input' in python 3.
    # Neither does 'input' in python2.
    if sys.platform.startswith('win'):
        raw_input('Press ENTER to exit:')
    # end if
    sys.exit(1)
# end if


from src.platform import platformwise_exit

# Firstly check if ve just need to print version or help message
if '-v' in sys.argv[1:] or '--version' in sys.argv[1:] or '-version' in sys.argv[1:]:
    print(__version__)
    platformwise_exit(0)
# end if


import src.print_help

if len(sys.argv) == 1 \
   or '-h' in sys.argv[1:] \
   or '--help' in sys.argv[1:] \
   or '-help' in sys.argv[1:]:
    src.print_help.print_help(__version__, __last_update_date__)
    platformwise_exit(0)
# end if


import os

print('\n  == {} v{} ==\n'.format(os.path.basename(__file__), __version__))


from src.main import main

if __name__ == '__main__':
    result_status = main()
    platformwise_exit(result_status)
# end if

