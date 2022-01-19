
import os
import argparse

from src.arguments import KromsatelArgs


def parse_args():
    argparse_args = _parse_command_line()
    kromsatel_args = KromsatelArgs(argparse_args)
    return kromsatel_args
# end def


def _parse_command_line():

    parser = argparse.ArgumentParser()

    # parser.add_argument(
    #     '-h',
    #     '--help',
    #     help='TODO',
    #     required=False,
    #     action="store_true"
    # )

    # parser.add_argument(
    #     '-v',
    #     '--version',
    #     help='TODO',
    #     required=False,
    #     action="store_true"
    # )

    parser.add_argument(
        '-1',
        '--reads-R1',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '-2',
        '--reads-R2',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '-l',
        '--reads-long',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '-p',
        '--primers',
        help='TODO',
        required=True
    )

    parser.add_argument(
        '-r',
        '--reference',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '-o',
        '--outdir',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '-s',
        '--split-output',
        help='TODO',
        required=False,
        action='store_true'
    )

    parser.add_argument(
        '-m',
        '--min-len',
        help='TODO',
        required=False,
        type=int
    )

    parser.add_argument(
        '-t',
        '--threads',
        help='TODO',
        required=False,
        type=int
    )

    parser.add_argument(
        '-c',
        '--chunk-size',
        help='TODO',
        required=False,
        type=int
    )

    parser.add_argument(
        '-k',
        '--blast-task',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '--crop-len',
        help='TODO',
        required=False
    )

    parser.add_argument(
        '--primer-5ext',
        help='TODO',
        required=False,
        type=int
    )

    parser.add_argument(
        '--use-index',
        help='TODO',
        required=False
    )

    args = parser.parse_args()

    return args
# end def
