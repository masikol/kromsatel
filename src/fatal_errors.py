
class FatalError(Exception):
    pass
# end class


class InvalidFastqError(FatalError):
    # For this error, we need to log errorneous line to a log file
    #   because this line may be too long for printing.
    # Hence these overridings and new fields.

    def __init__(self, /, msg_to_print, msg_to_log_only, *args, **kwargs):
        super().__init__(args, kwargs)
        self.msg_to_print = msg_to_print
        self.msg_to_log_only = msg_to_log_only
    # end def

    def __str__(self):
        return self.msg_to_print
    # end def

    def __repr__(self):
        return self.msg_to_print
    # end def
# end class
