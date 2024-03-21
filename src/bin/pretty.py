#!/usr/bin/env python3
"""Pretty outputs
"""

import argparse
import re

# Col scheme ------------------------------------------------------------------
COL_DNA = {
    "A": "red", "C": "gold", "G": "green",
    "T": "blue", "X": "grey", ".": "black"
}

COL_RES = {
    "A": "pink", "C": "gold", "D": "green",
    "d": "red", "E": "green", "e": "red",
    "F": "orange", "G": "magenta", "H": "red",
    "h": "red", "j": "red", "I": "pink",
    "K": "red", "L": "pink", "M": "pink",
    "N": "blue", "P": "magenta", "Q": "blue",
    "R": "red", "S": "blue", "T": "blue",
    "V": "pink", "W": "orange", "Y": "orange",
    "X": "grey", ".": "black"
}

# Functions -------------------------------------------------------------------

def clean_seq(full_seq, seq_prime, start, end):
    "return add gap for the missing positions"
    GAP = "."
    nb_gap = seq_prime.count("-")
    len_full = len(full_seq)
    right = len_full - end
    left = start
    return full_seq[:left]+seq_prime+full_seq[end:], left*GAP+seq_prime+right*GAP

def complement(seq):
    "find the complement"
    comp = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C'}
    return "".join([comp[el] for el in seq])


def color_seq(seq, colors, bold=True):
    "color if amino acid line"
    if bold:
        col = lambda char: '<b><font color="{}">{}</font></b>'.format(
            colors[char], char) if char in colors else char
    else:
        col = lambda char: '<font color="{}">{}</font>'.format(
            colors[char], char) if char in colors else char
    return "".join([col(char) for char in seq])

def pretty_print(frame_list, full_seqs, residues, res_prime, dna, html):
    "pretty print"
    by, by_dna = 26, 26 * 3
    names = ['X', 'Y']
    outputs = ""
    bold = lambda el: ("<b>" + el + "</b>" if html else el)
    color_res = (lambda s: color_seq(s, COL_RES, True)) if html else (lambda s: s)
    color_dna = (lambda s: color_seq(s, COL_DNA, False)) if html else (lambda s: s)
    color_res_seq = lambda s: "".join(map(color_res, s))
    color_dna_seq = lambda s: "".join(map(color_dna, s))

    outputs += "\n"
    for name_full, seq_full, seq_full_p in full_seqs:
        outputs += f">{name_full}" + "\n"
        outputs += color_res(seq_full) + "\n"
        outputs += color_res(seq_full_p) + "\n\n"
    outputs += "\n"

    upper_part = [(frame, name, seq_p[::-1], seq_i[::-1]) for frame, name, seq_p, seq_i in zip(frame_list, names, res_prime, residues) if frame < 0]
    down_part = [(frame, name, seq_p, seq_i) for frame, name, seq_p, seq_i in zip(frame_list, names, res_prime, residues) if frame > 0]

    for i in range(int(len(dna) / by_dna + 1)):
        for (frame, name, seq_p, seq_i) in upper_part:
            seq_i_ = seq_i[i * by:(i + 1) * by]
            seq_p_ = seq_p[i * by:(i + 1) * by]
            outputs += "{} ".format(name) + " " * abs(frame) + "  ".join(
                map(color_res, list(seq_i_))) + "\n"
            outputs += "{}'".format(name) + " " * abs(frame) + "  ".join(
                map(color_res, list(seq_p_))) + "\n"
        dna_i_ = dna[i * by_dna:(i + 1) * by_dna + 1]
        outputs += " " * 2 + "".join(
            map(color_dna,
                complement(dna_i_))) + "\n"
        outputs += " " * 2 + "".join(map(color_dna, dna_i_)) + "\n"
        for (frame, name, seq_p, seq_i) in down_part:
            seq_i_ = seq_i[i * by:(i + 1) * by]
            seq_p_ = seq_p[i * by:(i + 1) * by]
            outputs += "{}'".format(name) + " " * abs(frame) + "  ".join(
                map(color_res, list(seq_p_))) + "\n"
            outputs += "{} ".format(name) + " " * abs(frame) + "  ".join(
                map(color_res, list(seq_i_))) + "\n"
        outputs += "\n"
    return outputs


def to_html(output, out_file):
    with open(out_file, "w") as out:
        out.write("""
<html>
<pre>
<body>
{}
</body>
<pre>
</html>
        """.format(output))

def read_output(bin_out):
    "read the output from the binary"
    parm_reg = re.compile("# ([^ ]+) = ([^ ]+)")
    parm_type = {'NAMEX': str, 'NAMEY': str, 'GAP_PEN': float, 'GAP_OPEN': float, 'MAT': str, 'FRAME': int,'STARTX': int, 'STARTY': int, 'ENDX': int, 'ENDY': int, 'SCORE': float, 'DCA': str}
    parameters, sequences = {}, {}
    for l in bin_out:
        # SCORES and running infos
        if l.startswith("#"):
            parm_name, parm_val = parm_reg.search(l.strip()).groups()
            parameters[parm_name] = parm_type[parm_name](parm_val)
        elif not l == "":
            # read sequences
            if l.startswith(">"):
                name = l[1:].strip()
                sequences[name] = ""
            else:
                sequences[name] += l.strip()
    return parameters, sequences


def parse_arguments():
    "Parse command line"
    parser = argparse.ArgumentParser('pretty output')
    parser.add_argument("over_out", help="output from overgenes binaries")
    parser.add_argument('-o', '--out', help="output html")
    return parser.parse_args()


def main():
    args = parse_arguments()
    parameters, sequences = read_output(open(args.over_out))
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

    if args.out:
        output_print += pretty_print(frame_list, full_seq_l, residues, res_prime, dna, True)
        to_html(output_print, args.out)
    else:
        output_print += pretty_print(frame_list, full_seq_l, residues, res_prime, dna, False)
        print(output_print)

if __name__ == '__main__':
    main()
