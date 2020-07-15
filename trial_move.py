"""
Euler角形式で氷の構造を読みこみ、指定量のNH4とFをドープする。

結果はpickleとして出力。
"""

import sys
import random
from logging import getLogger, basicConfig, DEBUG, INFO
from math import sin, cos, pi
import numpy as np
import pairlist as pl
import pickle
from dope import trial_move, pick_anions, pick_cations



debug=True
if debug:
    basicConfig(level=DEBUG,
                        format="%(asctime)s %(levelname)s %(message)s")
else:
    basicConfig(level=INFO,
                        format="%(levelname)s %(message)s")

infile = sys.argv[1]
nMove  = int(sys.argv[2])
outfile = sys.argv[3]

with open(infile, "rb") as f:
    cellmat, rpos, anions, cations, G = pickle.load(f)
# anions, cationsは使わない。

nWater = rpos.shape[0]

while nMove > 0:
    result = trial_move(G)
    if result is not None:
        # accept the move
        G = result
        nMove -= 1

anions = pick_anions(G)
cations = pick_cations(G)
assert len(anions) == len(cations)

with open(outfile, "wb") as f:
    pickle.dump([cellmat, rpos, anions, cations, G], f)
