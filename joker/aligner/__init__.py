#!/usr/bin/env python3
# coding: utf-8

import weakref

from joker.aligner.algorithm import Aligner
from joker.aligner.utility import load_submat

__version__ = '0.0.3'

# cache a few aligners
_wref_aligners = weakref.WeakValueDictionary()


def create_aligner(submat='blosum60', rho=13, sigma=1):
    return Aligner.from_submat(*load_submat(submat), rho=rho, sigma=sigma)


def get_aligner(submat='blosum60', rho=13, sigma=1):
    key = submat, rho, sigma
    try:
        return _wref_aligners[key]
    except KeyError:
        pass
    a = create_aligner(submat=submat, rho=rho)
    _wref_aligners[key] = a
    return a
