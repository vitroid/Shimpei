"""
Euler角形式で氷の構造を読みこみ、指定量のNH4とFをドープする。

結果はpickleとして出力。
"""

import random
from logging import getLogger, basicConfig, DEBUG, INFO
import networkx as nx


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
    for i in G.successors(target):
        if G.in_degree(i) == 4:
            # if the neighbor is also an anion,
            logger.debug('Cancelled by conflict.')
            return False
    return True


def path2edges(path):
    return [(x,y) for x,y in zip(path, path[1:])]

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
    Embed nIon pairs of cations and anions in an ice represented by G
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



def move_anion(G, anion):
    logger = getLogger()
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
    return G


def move_cation(G, cation):
    logger = getLogger()
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
