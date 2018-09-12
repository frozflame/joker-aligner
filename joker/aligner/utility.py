#!/usr/bin/env python3
# coding: utf-8

from __future__ import unicode_literals, print_function

import os

import joker.aligner
import numpy as np


def locate_submat(name):
    d = os.path.split(joker.aligner.__file__)[0]
    possible_paths = [
        os.path.join(d, 'matrix', 'bioinfo', name),
        os.path.join(d, 'matrix', name),
        name,
    ]

    for p in possible_paths:
        if os.path.isfile(p):
            return p
    errmsg = 'cannot find substitue matrix "{}"'.format(name)
    raise IOError(errmsg)


def load_submat(path):
    if '/' not in path:
        path = locate_submat(path)

    with open(path) as infile:
        ichars = []
        jstring = None
        jsize = 0
        scores = []

        for line in infile:
            line = line.strip()
            if line and not line.startswith('#'):
                # remove whitespaces
                jstring = ''.join(line.split())
                jsize = len(jstring) + 1
                break

        if jstring is None:
            raise ValueError('bad substitution matrix')

        for line in infile:
            items = line.split()
            if not items:
                continue

            if len(items) != jsize:
                raise ValueError('bad substitution matrix')
            ichars.append(items[0])
            scores.append([int(s) for s in items[1:]])

    istring = ''.join(ichars)
    submatr = np.array(scores, dtype=int)
    return istring, jstring, submatr


