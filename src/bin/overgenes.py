#!/usr/bin/env python3
"""Run the overgene program binary:
"""

import re
import argparse
import subprocess
from os.path import dirname, abspath
from pretty import pretty_print, to_html, clean_seq, read_output


def parse_arguments():
    "Parse command line"
    parser = argparse.ArgumentParser('overlapping gene design or search')
    parser.add_argument("mode", help="local or global search (can be used for design)", choices=("local", "global"))
    parser.add_argument("fasta_x", help="input fasta files (one sequence per file)")
    parser.add_argument("fasta_y", help="input fasta files (one sequence per file)")
    parser.add_argument('-f', '--frame', help="overlapping reading frame", type=int, default=-2, choices=(-2, -1, 0, 1, 2, 3))
    parser.add_argument('-o', '--out', help="output file")
    parser.add_argument('-x', '--dca_x', help="dca score x")
    parser.add_argument('-y', '--dca_y', help="dca score x")
    parser.add_argument('-gp', '--gap_pen', help="gap penalty", default=-2.0, type=float)
    parser.add_argument('-go', '--gap_open', help="gap open", default=-16.0, type=float)
    parser.add_argument('-m', '--matrix', default="ident", choices=("ident", "blosum62", "blosum80", "blosum90", "vtml200"))
    parser.add_argument('--nopretty', help="use identity score", default=True, action="store_false")
    parser.add_argument('--raw', default=False, action="store_true", help="print sequences in fasta format")
    return parser.parse_args()


def main():
    args = parse_arguments()

    if args.mode.lower() == "local":
        binary = dirname(abspath(__file__)) + "/local"
    elif args.mode.lower() == "global":
        binary = dirname(abspath(__file__)) + "/global"

    matrix_type = args.matrix

    cmd_line = "{} {} {} --phase={} -s {} --gap_cost={} --gap_open={}".format(binary, args.fasta_x,
                                                                              args.fasta_y, args.frame,
                                                                              matrix_type,
                                                                              args.gap_pen,
                                                                              args.gap_open)
    if args.dca_x and args.dca_y:
        cmd_line += f" -x {args.dca_x} -y {args.dca_y}"

    output = subprocess.check_output(cmd_line, shell=True).decode('UTF-8')
    parameters, sequences = read_output(output.split("\n"))
    score = parameters['SCORE']
    seq_x, seq_y = sequences['SEQ_X_FULL'], sequences['SEQ_Y_FULL']
    seq_xi, seq_yi = sequences['SEQ_X_INIT'], sequences['SEQ_Y_INIT']
    seq_xp, seq_yp = sequences['SEQ_X_PRIME'], sequences['SEQ_Y_PRIME']
    dna = sequences['DNA']
    start_x, start_y = parameters['STARTX'], parameters['STARTY']
    end_x, end_y = parameters['ENDX'], parameters['ENDY']
    name_x, name_y = parameters['NAMEX'], parameters['NAMEY']
    frame = parameters['FRAME']
    frame_list = {0: [1, -1], 1: [1, 2], 2: [2, 1], 3: [1, 1], -1: [2, -1], -2: [1, -2]}[frame]

    if end_x < len(seq_x) and start_x > 0:
        clean_x_f, clean_xp = clean_seq(seq_x, seq_xp, start_x, end_x)
    else:
        clean_x_f, clean_xp = seq_xi, seq_xp

    if end_y < len(seq_y) and start_y > 0:
        clean_y_f, clean_yp = clean_seq(seq_y, seq_yp, start_y, end_y)
    else:
        clean_y_f, clean_yp = seq_yi, seq_yp

    name_x += "/{}-{}".format(start_x, end_x)
    name_y += "/{}-{}".format(start_y, end_y)

    residues = [list(seq_xi), list(seq_yi)]
    res_prime = [list(seq_xp), list(seq_yp)]
    full_seq_l = [(name_x, clean_x_f, clean_xp), (name_y, clean_y_f, clean_yp)]

    output_print = ""
    for parm, parm_v in parameters.items():
        output_print += f"# {parm} {parm_v}\n"

    if args.out is not None:
        output_print += pretty_print(frame_list, full_seq_l, residues, res_prime, dna, True)
        to_html(output_print, args.out)
    elif args.raw:
        output_print += f">{name_x}\n{clean_x_f}\n{clean_xp}\n"
        output_print += f">{name_y}\n{clean_y_f}\n{clean_yp}\n"
        output_print += f">DNA\n{dna}\n"
        print(output)
    else:
        output_print += pretty_print(frame_list, full_seq_l, residues, res_prime, dna, False)
        print(output_print)

if __name__ == '__main__':
    main()
