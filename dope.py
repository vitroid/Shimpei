"""
Euler角形式で氷の構造を読みこみ、指定量のNH4とFをドープする。

結果はpickleとして出力。
"""

import sys
import random
from logging import getLogger, basicConfig, DEBUG, INFO
from math import sin, cos, pi
import numpy as np
import networkx as nx
import pairlist as pl
import pickle

def load_nx3a(file):
    """
    重心位置とオイラー角を読みこむ。
    """
    coms = []
    while True:
        line = file.readline()
        if len(line) == 0:
            break
        if "@BOX3" in line:
            line = file.readline()
            cell = np.array([float(x) for x in line.split()])
        elif "@NX3A" in line:
            line = file.readline()
            nmol = int(line)
            for i in range(nmol):
                line = file.readline()
                comeu = [float(x) for x in line.split()[:6]]
                coms.append(comeu)
    return np.array(coms), cell

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

def hbn(rcom, cellmat, R, water, icetest=True):
    """
    分子の位置と配向から水素結合グラフを生成する
    """
    dg = nx.DiGraph()
    for i, j in pl.pairs_iter(rcom, 3.2, cellmat, distance=False):
        dij = rcom[j] - rcom[i]
        dij -= np.floor(dij + 0.5)
        dij = dij @ cellmat
        wi = water @ R[i]
        wj = water @ R[j] + dij
        dmin = 2.3**2 # ice 16のために長めにした。
        hb   = None
        for k in (0,1):
            dOH = wi[k] - wj[2]
            L2 = dOH @ dOH
            if L2 < dmin:
                dmin = L2
                hb   = (i, j) # i to j
            dOH = wj[k] - wi[2]
            L2 = dOH @ dOH
            if L2 < dmin:
                dmin = L2
                hb   = (j, i) # j to i
        if hb is not None:
            dg.add_edge(*hb)
    if icetest:
        for x in dg.out_degree(dg):
            assert x[1]==2, x
    return dg

def cationizable(G, target):
    """
    target位置にある水分子がcationに変換できるかどうかを判定する。
    """
    logger = getLogger()
    if G.out_degree(target) != 2:
        return False
    for i in G.predecessors(target):
        if G.out_degree(i) == 4:
            # if the neighbor is also a cation,
            logger.debug('Cancelled by conflict.')
            return False
    return True

def anionizable(G, target):
    logger = getLogger()
    if G.out_degree(target) != 2:
        return False
    for j in G.successors(target):
        if G.in_degree(i) == 4:
            # if the neighbor is also an anion,
            logger.debug('Cancelled by conflict.')
            return False
    return True

def quat2rotmat(q):
    """
    Quaternionを回転行列(後置形式)にする。
    """
    a, b, c, d = q
    sp11 = (a * a + b * b - (c * c + d * d))
    sp12 = -2.0 * (a * d + b * c)
    sp13 = 2.0 * (b * d - a * c)
    sp21 = 2.0 * (a * d - b * c)
    sp22 = a * a + c * c - (b * b + d * d)
    sp23 = -2.0 * (a * b + c * d)
    sp31 = 2.0 * (a * c + b * d)
    sp32 = 2.0 * (a * b - c * d)
    sp33 = a * a + d * d - (b * b + c * c)
    return np.array([[sp11, sp12, sp13], [sp21, sp22, sp23], [sp31, sp32, sp33]]).T

def euler2quat(e):
    """
    Euler角をQuaternionにする。
    """
    ea, eb, ec = e
    a = cos(ea / 2) * cos((ec + eb) / 2)
    b = sin(ea / 2) * cos((ec - eb) / 2)
    c = sin(ea / 2) * sin((ec - eb) / 2)
    d = cos(ea / 2) * sin((ec + eb) / 2)
    return np.array((a, b, c, d))

def path2edges(path):
    return [tuple(path[i:i+2]) for i in range(len(path)-1)]

def remove_path(G, path):
    """
    Gからpathを消す。消したパスをDiGraphとして返す。
    """
    rG = G.edge_subgraph(edges)

    for i in range(len(path)-1):
        a,b = path[i:i+2]

def invert_edge(G, from_, to_):
    """
    invert an edge
    """
    assert G.has_edge(from_, to_)
    G.remove_edge(from_, to_)
    G.add_edge(to_, from_)

def bulkdope(G, nIon):
    """
    水の水素結合ネットワークGの中に、nIon対のイオンを埋めこむ。
    それにより、水素結合ネットワークのトポロジーが変化する。
    """
    logger = getLogger()
    nWater = G.number_of_nodes()
    logger.info(nWater)
    anions = set()
    cations = set()
    N = nIon
    while N > 0:
        # アニオンを置く場所をランダムに選ぶ。
        while True:
            anion = random.randint(0,nWater-1)
            if anionizable(G, anion):
                break
        # カチオンを置く場所をランダムに選ぶ。
        while True:
            cation = random.randint(0,nWater-1)
            if cation != anion and cationizable(G, cation):
                break
        # 水2つをイオン2つに置き換えたいなら、2本のパスを反転させる必要がある。
        try:
            path1 = nx.shortest_path(G, anion, cation)
        except nx.exception.NetworkXNoPath:
            # 1つもパスがないのならあきらめる。
            continue
        logger.debug("Path1: {0}".format(path1))
        # homodromic path上にイオンが含まれる心配はない。
        # pathをedgeの集合に変換
        edges1 = path2edges(path1)
        # 2本目の経路は、1本目の経路を除去してから探す。
        # パスの部分だけ抽出したグラフ
        rG1    = G.edge_subgraph(edges1).copy()
        # パスを除去したグラフ
        xG = G.copy()
        xG.remove_edges_from(edges1)
        try:
            path2 = nx.shortest_path(xG, anion, cation)
        except nx.exception.NetworkXNoPath:
            # パスがないのならあきらめる。
            continue
        logger.debug("Path2: {0}".format(path2))
        edges2 = path2edges(path2)
        for f,t in edges1:
            invert_edge(G, f, t)
        for f,t in edges2:
            invert_edge(G, f, t)
        # イオンの周囲の水素結合の向きを確認
        assert G.in_degree(anion) == 4
        assert G.out_degree(cation) == 4
        anions.add(anion)
        cations.add(cation)
        N -= 1
        logger.info(N)
    return G, anions, cations


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
# 2番目の引数はパーセンテージ。10を指定するとanion, cationそれぞれ10%、水が80%
percent = float(sys.argv[2])
# 出力ファイルは今のところバイナリ。
outfile = sys.argv[3]

# 氷の構造を読みこむ
waters, cell = load_nx3a(open(infile))
pos   = waters[:,:3]
euler = waters[:,3:6]
# セル行列 (行ごとにa,b,c)
cellmat = np.diag(cell)
# セルの逆行列
celli   = np.linalg.inv(cellmat)
# 水分子の総数
nWater = pos.shape[0]
# セル内相対座標
rpos = pos @ celli
# 水分子の配向行列
R = np.zeros([nWater,3,3])
for i in range(nWater):
    R[i] = quat2rotmat(euler2quat(euler[i]))
# 水分子の原子位置
water = tip4picesites()
# 水素結合ネットワーク(有向グラフ)
G = hbn(rpos, cellmat, R, water)

# イオンの割合 10だとアニオン10% カチオンも10%、水は80%
nIon = int(nWater * percent / 100)

# ドープする。
G, anions, cations = bulkdope(G, nIon)

# できたので保存する。
# 必要な情報は、anionとcationの位置、格子点の位置、水のグラフ
with open(outfile, "wb") as f:
    pickle.dump([cellmat, rpos, anions, cations, G], f)
