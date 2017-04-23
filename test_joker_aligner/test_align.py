#!/usr/bin/env python3
# coding: utf-8
from __future__ import unicode_literals, print_function

import itertools
import sys
# import six
import time

from joker.aligner.calculate import Aligner
from joker.aligner.submatrix import load_submat

# from molbiox.frame.environ import locate_tests
# from molbiox.frame.testing import Timer
# from molbiox.io import submat, fasta


# sp|P37744|RMLA1_ECOLI Glucose-1-phosphate thymidylyltransferase 1
# OS=Escherichia coli (strain K12) GN=rfbA PE=1 SV=2
iseq = """
    MKMRKGIILAGGSGTRLYPVTMAVSKQLLPIYDKPMIYYPLSTLMLAGIRDILIISTPQD
    TPRFQQLLGDGSQWGLNLQYKVQPSPDGLAQAFIIGEEFIGGDDCALVLGDNIFYGHDLP
    KLMEAAVNKESGATVFAYHVNDPERYGVVEFDKNGTAISLEEKPLEPKSNYAVTGLYFYD
    NDVVQMAKNLKPSARGELEITDINRIYLEQGRLSVAMMGRGYAWLDTGTHQSLIEASNFI
    ATIEERQGLKVSCPEEIAFRKGFIDVEQVRKLAVPLIKNNYGQYLYKMTKDSN
    """


# sp|P55255|RMLA_NEIMB Glucose-1-phosphate thymidylyltransferase
# OS=Neisseria meningitidis serogroup B (strain MC58) GN=rmlA1 PE=3 SV=2
jseq = """
    MKGIILAGGSGTRLYPITRGVSKQLLPVYDKPMIYYPLSVLMLAGIRDILVITAPEDNAS
    FKRLLGDGSDFGISISYAVQPSPDGLAQAFIIGEEFIGNDNVCLVLGDNIFYGQSFTQTL
    KQAAAQTHGATVFAYQVKNPERFGVVEFNENFRAVSIEEKPQRPKSDWAVTGLYFYDNRA
    VEFAKQLKPSARGELEITDLNRMYLEDGSLSVQILGRGFAWLDTGTHESLHEAASFVQTV
    QNIQNLHIACLEEIAWRNGWLSDEKLEELARPMAKNQYGQYLLRLLKK
    """


iseq = ''.join(iseq.split())
jseq = ''.join(jseq.split())


score_type = 'int64'
repeat = 1000


def chunkwize_parallel(chunksize, *args):
    # args are strings or lists
    chunksize = int(chunksize)
    for i in itertools.count(0):
        r = [s[i*chunksize:(i+1)*chunksize] for s in args]
        if any(r):
            yield r
        else:
            raise StopIteration


def test_align():
    istring, jstring, submatr = load_submat('pam200')
    aligner = Aligner.from_submatr(istring, jstring, submatr)

    tm = time.time()
    for i in range(repeat):
        matrx, istring, jstring = aligner.align(iseq, jseq, backtrack=True)
    tm = time.time() - tm
    msg = '{} loops, {:.4f} sec per loop'.format(repeat, tm / repeat)
    print(msg)

    mstring = Aligner.make_match_string(istring, jstring)

    for items in chunkwize_parallel(60, istring, mstring, jstring):
        bunch = b'\n'.join(items).decode('ascii')
        print(bunch, end='\n\n', file=sys.stderr)


if __name__ == '__main__':
    test_align()



