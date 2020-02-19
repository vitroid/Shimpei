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
import yaplotlib as yp
import itertools as it

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


def draw_water(water, povray=False):
    s = ""
    if povray:
        s += pov.Sphere(water[0], r="RH", material="MATH")
        s += pov.Sphere(water[1], r="RH", material="MATH")
        s += pov.Sphere(water[2], r="RO", material="MATO")
        s += pov.Cylinder(water[0], water[2], r="ROH", material="MATOH")
        s += pov.Cylinder(water[1], water[2], r="ROH", material="MATOH")
    else:
        s += yp.Layer(1)
        s += yp.Color(5)
        s += yp.Size(0.2)
        s += yp.Circle(water[0])
        s += yp.Circle(water[1])
        s += yp.Color(4)
        s += yp.Size(0.4)
        s += yp.Circle(water[2])
        s += yp.Color(2)
        d0 = water[0] - water[2]
        L0 = np.linalg.norm(d0)
        s += yp.Line(water[2] + d0/L0*0.4, water[0] - d0/L0*0.2)
        d1 = water[1] - water[2]
        L1 = np.linalg.norm(d1)
        s += yp.Line(water[2] + d1/L1*0.4, water[1] - d1/L1*0.2)
        # draw a tetrahedron
        y = water[1] - water[0]
        z = (water[1] + water[0]) / 2 - water[2]
        x = np.cross(y,z)
        x /= np.linalg.norm(x)
        y /= np.linalg.norm(y)
        z /= np.linalg.norm(z)
        com = (water[1] + water[0] + water[2]*16)/18
        R = 2.76/2
        a = y*(2/3)**0.5 + z*(1/3)**0.5
        b = -y*(2/3)**0.5 + z*(1/3)**0.5
        c = x*(2/3)**0.5 - z*(1/3)**0.5
        d = -x*(2/3)**0.5 - z*(1/3)**0.5
        s += yp.Layer(2)
        for e,f in it.combinations((a,b,c,d), 2):
            s += yp.Line(com+e*R,com+f*R)
    return s

def to_yaplot(cellmat, rpos, anions, cations, G, cation="N", anion="F"):
    """
    cellmat: セル行列(後置記法)
    rpos:    分子のセル内相対位置
    anions:  アニオン番号
    cations: カチオン番号
    G:       水素結合ネットワーク
    """
    # 水分子の原子位置
    water = tip4picesites()

    s = ""
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
            s += draw_water(intra)
    s += yp.Layer(3)
    s += yp.Size(0.5)
    s += yp.Color(6)
    for mol in G:
        if mol in anions:
            pos = rpos[mol] @ cellmat
            s += yp.Circle(pos)
    s += yp.Size(0.3)
    s += yp.Color(7)
    for mol in G:
        if mol in cations:
            pos = rpos[mol] @ cellmat
            s += yp.Circle(pos)
    s += yp.NewPage()
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

infile = sys.argv[1]
with open(infile, "rb") as f:
    cellmat, rpos, anions, cations, G = pickle.load(f)


s = to_yaplot(cellmat, rpos, anions, cations, G)
print(s)
