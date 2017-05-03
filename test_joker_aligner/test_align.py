#!/usr/bin/env python3
# coding: utf-8
from __future__ import unicode_literals, print_function

import itertools
import sys
import time

from joker.aligner import align_simple

try:
    from joker.textmanip.iterative import chunkwize_parallel
except ImportError:
    def chunkwize_parallel(chunksize, *args):
        # args are strings or lists
        chunksize = int(chunksize)
        for i in itertools.count(0):
            r = [s[i*chunksize:(i+1)*chunksize] for s in args]
            if any(r):
                yield r
            else:
                raise StopIteration


# sp|P37744|RMLA1_ECOLI Glucose-1-phosphate thymidylyltransferase 1
# OS=Escherichia coli (strain K12) GN=rfbA PE=1 SV=2
i_sequence = """
    MKMRKGIILAGGSGTRLYPVTMAVSKQLLPIYDKPMIYYPLSTLMLAGIRDILIISTPQD
    TPRFQQLLGDGSQWGLNLQYKVQPSPDGLAQAFIIGEEFIGGDDCALVLGDNIFYGHDLP
    KLMEAAVNKESGATVFAYHVNDPERYGVVEFDKNGTAISLEEKPLEPKSNYAVTGLYFYD
    NDVVQMAKNLKPSARGELEITDINRIYLEQGRLSVAMMGRGYAWLDTGTHQSLIEASNFI
    ATIEERQGLKVSCPEEIAFRKGFIDVEQVRKLAVPLIKNNYGQYLYKMTKDSN
    """

# iseq = iseq[60:120]


# sp|P55255|RMLA_NEIMB Glucose-1-phosphate thymidylyltransferase
# OS=Neisseria meningitidis serogroup B (strain MC58) GN=rmlA1 PE=3 SV=2
j_sequence = """
    MKGIILAGGSGTRLYPITRGVSKQLLPVYDKPMIYYPLSVLMLAGIRDILVITAPEDNAS
    FKRLLGDGSDFGISISYAVQPSPDGLAQAFIIGEEFIGNDNVCLVLGDNIFYGQSFTQTL
    KQAAAQTHGATVFAYQVKNPERFGVVEFNENFRAVSIEEKPQRPKSDWAVTGLYFYDNRA
    VEFAKQLKPSARGELEITDLNRMYLEDGSLSVQILGRGFAWLDTGTHESLHEAASFVQTV
    QNIQNLHIACLEEIAWRNGWLSDEKLEELARPMAKNQYGQYLLRLLKK
    """


i_sequence = ''.join(i_sequence.split())
j_sequence = ''.join(j_sequence.split())


def demo(iseq, jseq, quiet=False):
    score, iseq, jseq, mseq = align_simple(iseq, jseq)
    for items in chunkwize_parallel(60, iseq, mseq, jseq):
        bunch = '\n'.join(items)
        if not quiet:
            print(bunch, end='\n\n', file=sys.stderr)
    if not quiet:
        print('score:', score)
    return score, iseq, jseq, mseq


def benchmark(loop=100):
    tm = time.time()
    for i in range(loop):
        _ = align_simple(i_sequence, j_sequence)
    tm = time.time() - tm
    tm_perloop = 1000 * tm / loop
    msg = '{} loops, avg {:.2f} millisec per loop'.format(loop, tm_perloop)
    print(msg)


def test_align():
    # test_align()
    demo(i_sequence, j_sequence)
    benchmark()


if __name__ == '__main__':
    test_align()


