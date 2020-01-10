#!/usr/bin/env python3
# coding: utf-8

from __future__ import division, print_function

# from joker.aligner import visualize
# from joker.aligner.compute import Aligner
from joker.aligner.utils import load_submat


def test_load_submat():
    i, j, s = load_submat('pam200')
    print(i)
    print(j)
    print(s)


def test_load_alnum1():
    i, j, s = load_submat('ENSUB')
    print(i)
    print(j)
    print(s)


if __name__ == '__main__':
    test_load_submat()
    test_load_alnum1()
