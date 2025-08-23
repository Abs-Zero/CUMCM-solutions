from numpy import pi, sin, cos, log, sqrt, arctan2
import numpy as np
from numba import jit
from matplotlib import pyplot as plt
import matplotlib as mpl
from openpyxl.styles import Alignment, Font
import openpyxl as px


mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False
plt.style.use("ggplot")


global V_HEAD, N_QUEUE, L0, L1, R, D, B, TH_MIN, S_MIN, BD2
V_HEAD = 1
N_QUEUE = 223
L0, L1 = 2.86, 1.65
R = 4.5
D = 1.7
B = D / 2 / pi
TH_MIN = R / B
S_MIN = B / 2 * (TH_MIN * sqrt(1 + TH_MIN**2) + log(TH_MIN + sqrt(1 + TH_MIN**2)))
BD2 = B / 2

global O1_X, O1_Y, O2_X, O2_Y, R1, R2, C1, C12, ALPHA_11, ALPHA_21
# 已知盘入螺线和盘出螺线与掉头区域交界处的极角 TH_MIN
# x,y 为两切点坐标，nx,ny 为两切点切线方向
# O1X,O1Y 为圆弧圆心
px, py = B * TH_MIN * cos(TH_MIN), B * TH_MIN * sin(TH_MIN)
nx, ny = -B * cos(TH_MIN) + TH_MIN * B * sin(TH_MIN), -B * sin(
    TH_MIN
) - TH_MIN * B * cos(TH_MIN)
nq = sqrt(nx * nx + ny * ny)
nx, ny = nx / nq, ny / nq
k1, k2 = 3, 2
r = (px * px + py * py) / ((k1 + k2) * (nx * py - ny * px))
R1, R2 = k1 * r, k2 * r
O1_X, O1_Y = px + ny * R1, py - nx * R1
O2_X, O2_Y = -px - ny * R2, -py + nx * R2
mx, my = (k2 * O1_X + k1 * O2_X) / (k1 + k2), (k2 * O1_Y + k1 * O2_Y) / (k1 + k2)

ALPHA_11 = arctan2(py - O1_Y, px - O1_X)
ALPHA_12 = arctan2(my - O1_Y, mx - O1_X)
ALPHA_21 = arctan2(O1_Y - my, O1_X - mx)
ALPHA_22 = arctan2(-py - O2_Y, -px - O2_X)

_BETA_1 = (ALPHA_11 - ALPHA_12) % (2 * pi)
_BETA_2 = (ALPHA_22 - ALPHA_21) % (2 * pi)

ALPHA_12 = ALPHA_11 - _BETA_1
ALPHA_22 = ALPHA_21 + _BETA_2

C1, C2 = R1 * _BETA_1, R2 * _BETA_2
C12 = C1 + C2


def plot_two_circle():
    linewidth = 2
    def plot_circle(x, y, r, s=0, e=2 * pi):
        theta = np.linspace(s, e, 100)
        x1 = r * cos(theta) + x
        x2 = r * sin(theta) + y
        plt.plot(x1, x2, color="black", linewidth=linewidth)

    def plot_single_segment(A, B):
        x1, y1 = A
        x2, y2 = B
        plt.plot(
            [x1, x2],
            [y1, y2],  # x坐标列表和y坐标列表
            color="black",  # 线段颜色
            linewidth=linewidth,  # 线宽
            linestyle="-",  # 线型（实线）
            marker="o",  # 端点标记（圆圈）
            markersize=6,
        )  # 标记大小

    def plot_vt_segment(A, B):
        n = B - A
        I = [[0, 1], [-1, 0]]
        n = I @ n
        C = B - n * 1
        D = B + n * 1
        plot_single_segment(C, D)
        return C, D

    def plot_points(**Ps):
        for k, v in Ps.items():
            plt.scatter(v[0], v[1], color="black", label=k)
            plt.text(
                v[0] + 0.3,
                v[1],
                k,
                color="black",
                fontsize=12,
            )

    def foot_of_perpendicular(P, A, B):
            px, py = P
            ax, ay = A
            bx, by = B
            
            # 向量AB和AP
            abx = bx - ax
            aby = by - ay
            apx = px - ax
            apy = py - ay
            
            # 计算投影参数t
            ab_dot_ab = abx **2 + aby** 2  # |AB|²
            if ab_dot_ab == 0:  # A和B重合时，垂足为A
                return (ax, ay)
            
            t = (apx * abx + apy * aby) / ab_dot_ab  # 点积除以|AB|²
            t_clamped = max(0, min(1, t))  # 约束t在[0,1]
            
            # 计算垂足坐标
            qx = ax + t_clamped * abx
            qy = ay + t_clamped * aby
            
            return np.array([qx, qy])

    O1 = np.array([O1_X, O1_Y])
    O2 = np.array([O2_X, O2_Y])
    M = (2 * O1 + 3 * O2) / 5
    A = O1 + np.array([R1 * cos(ALPHA_11), R1 * sin(ALPHA_11)])
    B = O2 + np.array([R2 * cos(ALPHA_22), R2 * sin(ALPHA_22)])
    
    
    plot_circle(O1[0], O1[1], R1, ALPHA_11, ALPHA_11 - C1 / R1)
    plot_circle(O2[0], O2[1], R2, ALPHA_21, ALPHA_21 + (C12 - C1) / R2)
    plot_single_segment(O1, M)
    plot_single_segment(O2, M)
    plot_single_segment(O1, A)
    C, D = plot_vt_segment(O1, A)
    H = foot_of_perpendicular(O2, C, D)
    plot_single_segment(H, O2)
    plot_single_segment(O2, B)
    plot_vt_segment(O2, B)
    plot_points(O1=O1, O2=O2, M=M, A=A, B=B, H=H)
    plt.axis("equal")
    plt.grid(False)
    plt.savefig(f"p4-示意图.svg")
    plt.clf()


plot_two_circle()
