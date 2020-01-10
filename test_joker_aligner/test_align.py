#!/usr/bin/env python3
# coding: utf-8
from __future__ import unicode_literals, print_function

import sys
import time

from joker.aligner import get_aligner
from joker.aligner.utils import chunkwize_parallel

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


def demostrate(iseq, jseq):
    aligner = get_aligner()
    alignment = aligner(iseq, jseq, backtrack=True)
    match_strings = (
        alignment.get_istr(),
        alignment.get_xstr(),
        alignment.get_jstr(),
    )
    for items in chunkwize_parallel(60, *match_strings):
        bunch = '\n'.join(items)
        print(bunch, end='\n\n', file=sys.stderr)
    print('score:', alignment.score)


def benchmark(loop=100):
    tm = time.time()
    aligner = get_aligner()
    for i in range(loop):
        _ = aligner(i_sequence, j_sequence)
    tm = time.time() - tm
    tm_perloop = 1000 * tm / loop
    msg = '{} loops, avg {:.2f} millisec per loop'.format(loop, tm_perloop)
    print(msg)


def test_align():
    # test_align()
    demostrate(i_sequence, j_sequence)
    benchmark()


def test():
    from joker.aligner import get_aligner, create_aligner
    s1 = 'KQAAAQTHGATVFAYQVKNPERFGVVEFNENFRAVS'
    s2 = 'KQAAAQTHGATVFGVVEFNENFzRAVSIEEKPQRPK'
    aligner = create_aligner('blosum60')
    # ali = aligner(s1, s2, backtrack=False)
    ali = aligner(s1, s2, backtrack=True)
    print(ali.get_istr())
    print(ali.get_xstr())
    print(ali.get_jstr())


if __name__ == '__main__':
    test_align()
    # test()
