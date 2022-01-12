
import sys


def platf_depend_exit(exit_code):
    if sys.platform.startswith('win'):
        input('Press ENTER to exit:')
    # end if
    sys.exit(exit_code)
# end def
