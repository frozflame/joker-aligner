#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

import sys

import numpy as np


def debug_trace(matrx):
    symbols = ord('-'), ord('|'), ord('\\'), ord('3')
    symbols = np.array(symbols, dtype='uint8')
    idx = matrx[:, :, 3]
    darr = symbols[idx]
    for row in darr:
        line = row.tostring().decode('ascii')
        print(line, file=sys.stderr)
