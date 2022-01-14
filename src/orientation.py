
class Orientation():

    LEFT  = 0
    RIGHT = 1

    _INVERT_DATA_STRUCT = (RIGHT, LEFT)

    @classmethod
    def invert(cls, orientation):
        return cls._INVERT_DATA_STRUCT[orientation]
    # end def
# end class


def get_read_orientation(alignment):

    orientation_left = alignment.align_strand_plus

    if orientation_left:
        return Orientation.LEFT
    else:
        return Orientation.RIGHT
    # end if
# end def
