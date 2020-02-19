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
    fix = G[from_][to_]['fixed']
    G.remove_edge(from_, to_)
    G.add_edge(to_, from_, fixed=fix)



def move_anion(G, anion):
    # 隣の水分子をさがす
    neis = list(G.predecessors(anion))
    random.shuffle(neis)
    assert len(neis) == 4
    found = False
    for water in neis:
        if G.in_degree(water) == 2:
            found = True
            break
    # となりに水分子がなければあきらめる
    if not found:
        return None

    # となりとの結合を除去したグラフ
    xG = G.copy()
    xG.remove_edge(water, anion)
    try:
        path = nx.shortest_path(xG, water, anion)
    except nx.exception.NetworkXNoPath:
        # パスがないのならあきらめる。
        return None
    logger.debug("Path: {0}".format(path))
    edges = path2edges(path)
    for f,t in edges:
        invert_edge(G, f, t)
    invert_edge(G, water, anion)
    # ややこしいので、ラベルを交換
    water, anion = anion, water
    # イオンと水の周囲の水素結合の向きを確認
    assert G.in_degree(water) == 2
    assert G.in_degree(anion) == 4
    # anionの周囲のHBは無条件で固定
    for nei in G.predecessors(anion):
        G[nei][anion]['fixed'] = True
    # 水の周囲については、
    # 水はもともとはanionだったので、結合はすべてfixedであった。
    # 2本の結合が反転し、outboundに変わった。
    # そのうち1つは新anionとの結合で無条件にfixedされる。
    # もう一方は、必ず水分子なので、fixedしない。
    for nei in G.successors(water):
        if nei != anion:
            G[water][nei]['fixed'] = False
    return G


def move_cation(G, cation):
    # 隣の水分子をさがす
    neis = list(G.successors(cation))
    random.shuffle(neis)
    assert len(neis) == 4
    found = False
    for water in neis:
        if G.in_degree(water) == 2:
            found = True
            break
    # となりに水分子がなければあきらめる
    if not found:
        return None

    # となりとの結合を除去したグラフ
    xG = G.copy()
    xG.remove_edge(cation, water)
    try:
        path = nx.shortest_path(xG, cation, water)
    except nx.exception.NetworkXNoPath:
        # パスがないのならあきらめる。
        return None
    logger.debug("Path: {0}".format(path))
    edges = path2edges(path)
    for f,t in edges:
        invert_edge(G, f, t)
    invert_edge(G, cation, water)
    # ややこしいので、ラベルを交換
    water, cation = cation, water
    # イオンと水の周囲の水素結合の向きを確認
    assert G.in_degree(water) == 2
    assert G.out_degree(cation) == 4
    # anionの周囲のHBは無条件で固定
    for nei in G.predecessors(cation):
        G[cation][nei]['fixed'] = True
    # 水の周囲については、
    # 水はもともとはcationだったので、結合はすべてfixedであった。
    # 2本の結合が反転し、inboundに変わった。
    # そのうち1つは新cationとの結合で無条件にfixedされる。
    # もう一方は、必ず水分子なので、fixedしない。
    for nei in G.predecessors(water):
        if nei != cation:
            G[nei][water]['fixed'] = False
    return G

def pick_anions(G):
    anions = []
    for i in G:
        if G.in_degree(i) == 4:
            anions.append(i)
    return anions

def pick_cations(G):
    cations = []
    for i in G:
        if G.out_degree(i) == 4:
            cations.append(i)
    return cations



def trial_move(G):
    """
    Gのなかのイオンを1つ選び、隣の格子点に移動する。失敗したらNoneを返す。
    """
    logger = getLogger()
    nWater = G.number_of_nodes()
    # イオンをランダムに選ぶ。
    while True:
        ion = random.randint(0,nWater-1)
        if G.in_degree(ion) == 4:
            # it is anion!
            return move_anion(G, ion)
        elif G.out_degree(ion) == 4:
            # it is cation!
            return move_cation(G, ion)


debug=True
if debug:
    basicConfig(level=DEBUG,
                        format="%(asctime)s %(levelname)s %(message)s")
else:
    basicConfig(level=INFO,
                        format="%(levelname)s %(message)s")
logger = getLogger()

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

with open(outfile, "wb") as f:
    pickle.dump([cellmat, rpos, anions, cations, G], f)
