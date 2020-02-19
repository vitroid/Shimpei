"""
Euler角形式で氷の構造を読みこみ、指定量のNH4とFをドープする。

結果はpickleとして出力。
"""

import sys
from logging import getLogger, basicConfig, DEBUG, INFO
import numpy as np
import networkx as nx
import pickle
from math import sin, cos, pi

def tip4picesites():
    """
    TIP4P/Ice Geometry
    """
    L1 = 0.9572
    L2 = 0.1577
    theta=104.52 * pi/180


    hy = L1*sin(theta/2)
    hz = L1*cos(theta/2)
    mz = L2

    # HHOM
    sites = np.array([[0.0, hy,  hz],
                      [0.0,-hy,  hz],
                      [0.0, 0.0, 0.0],
                      [0.0, 0.0, mz]])
    sites -= (sites[0]+sites[1]+sites[3]*0)/18
    return sites

def to_mdview(cellmat, rpos, anions, cations, G, cation="N", anion="F"):
    """
    cellmat: セル行列(後置記法)
    rpos:    分子のセル内相対位置
    anions:  アニオン番号
    cations: カチオン番号
    G:       水素結合ネットワーク
    """
    # 水分子の原子位置
    water = tip4picesites()
    AU    = 0.528

    atoms = []
    for mol in G:
        if mol not in anions and mol not in cations:
            # genuine water molecule.
            rO = rpos[mol]
            H1, H2 =  G[mol]
            rH1 = rpos[H1] - rO
            rH1 -= np.floor(rH1+0.5)
            rH2 = rpos[H2] - rO
            rH2 -= np.floor(rH2+0.5)
            # 実座標 (angstrom)
            pO = rO @ cellmat
            pH1 = rH1 @ cellmat
            pH2 = rH2 @ cellmat
            ey = pH1 - pH2
            ez = pH1 + pH2
            ex = np.cross(ey, ez)
            # 分子内座標方向の単位ベクトル
            ex /= np.linalg.norm(ex)
            ey /= np.linalg.norm(ey)
            ez /= np.linalg.norm(ez)
            # 回転行列
            R = np.array([ex,ey,ez])
            # 水分子の原子位置
            intra = water @ R + pO
            atoms.append(["O", intra[2]])
            atoms.append(["H", intra[0]])
            atoms.append(["H", intra[1]])
    for mol in G:
        if mol in anions:
            pos = rpos[mol] @ cellmat
            atoms.append([anion, pos])
    for mol in G:
        if mol in cations:
            pos = rpos[mol] @ cellmat
            atoms.append([cation, pos])
    s = "{0}\n".format(len(atoms))
    for name, pos in atoms:
        s += "{0} {1:.2f} {2:.2f} {3:.2f}\n".format(name, *(pos/AU)) # in atomic units
    return s

# debug messageの準備
debug=False
if debug:
    basicConfig(level=DEBUG,
                        format="%(asctime)s %(levelname)s %(message)s")
else:
    basicConfig(level=INFO,
                        format="%(levelname)s %(message)s")
logger = getLogger()

# 今のところ座標ファイルを読みこむことにしているが、
# ドーピングに必要なのは水素結合ネットワークのトポロジーのみ。
infile = sys.argv[1]
with open(infile, "rb") as f:
    cellmat, rpos, anions, cations, G = pickle.load(f)


s = to_mdview(cellmat, rpos, anions, cations, G)
print(s)
