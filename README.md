joker-aligner
=============

Global, local and overlap sequence alignment.

```python
from joker.aligner import get_aligner

s1 = 'MKMRGIILAGGSGTRLYPVTMAVSKQLLPIYDKPMIYYPLSTLMLAGIRDILIISTPQD'
s2 = 'MKGIKILAGGSGSSTRLYPITRGVSKQLLPVYDKPMLAGIRDILVITAPENASTTTTTT'

aligner = get_aligner('blosum62')
ali = aligner(s1, s2)



```
