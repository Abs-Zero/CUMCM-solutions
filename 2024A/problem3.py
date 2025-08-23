from numpy import pi, sin, cos, log, sqrt
import numpy as np
from numba import jit
from matplotlib import pyplot as plt
import matplotlib as mpl


mpl.rcParams["font.sans-serif"] = ["SimHei"]
mpl.rcParams["axes.unicode_minus"] = False
plt.style.use("ggplot")

N_QUEUE = 223  # 龙队节数
L0, L1 = 2.86, 1.65  # 龙头与龙身长度
H1, H2 = 0.275, 0.15  # 把手到边缘的距离

TH_MAX = 16 * 2 * pi  # 龙头初始极角
V_HEAD = 1  # 龙头速度
B_PRE = 1 / 2 / pi
S_MAX_PRE = 1 / 2 * (TH_MAX * sqrt(1 + TH_MAX**2) + log(TH_MAX + sqrt(1 + TH_MAX**2)))

R = 4.5
TH_MIN_PRE = R / B_PRE


@jit
def s_to_xy(s, b, s_max):
    def inver_s(s, b):
        # 螺线长的反函数，由自然坐标求极角
        s *= 2 / b
        th = sqrt(s)
        while True:
            th_new = th / 2 + (s - log(th + sqrt(1 + th**2))) / (2 * sqrt(1 + th**2))
            if th - th_new < 1e-12:
                break
            th = th_new
        return th_new

    s = s_max - s
    th = inver_s(s, b)
    return (
        b * th * cos(th),
        b * th * sin(th),
        -(cos(th) - th * sin(th)) / sqrt(1 + th**2),
        -(sin(th) + th * cos(th)) / sqrt(1 + th**2),
    )


@jit
def next_point(s_prev, x_prev, y_prev, length, b, s_max):
    # 返回下一个点的自然坐标，x坐标，y坐标，切向x分量，切向y分量
    s = s_prev - length
    while True:
        x, y, tx, ty = s_to_xy(s, b, s_max)
        d2 = (x - x_prev) ** 2 + (y - y_prev) ** 2
        L = d2 - length * sqrt(d2)
        Lp = 2 * ((x - x_prev) * tx + (y - y_prev) * ty)
        s_new = s - L / Lp

        if s - s_new < 1e-12:
            break
        s = s_new
    return s_new, x, y, tx, ty


@jit
def queue(r, d):
    # 返回龙头极径为r时，龙队各点信息
    b = d * B_PRE
    s_max = b * S_MAX_PRE
    th = r / b
    s = s_max - b / 2 * (th * (th + 1) + log(th + sqrt(1 + th**2)))
    ss, xs, ys = [0.0] * (N_QUEUE + 1), [0.0] * (N_QUEUE + 1), [0.0] * (N_QUEUE + 1)
    txs, tys = [0.0] * (N_QUEUE + 1), [0.0] * (N_QUEUE + 1)
    ss[0] = s
    xs[0], ys[0], txs[0], tys[0] = s_to_xy(ss[0], b, s_max)
    ss[1], xs[1], ys[1], txs[1], tys[1] = next_point(ss[0], xs[0], ys[0], L0, b, s_max)
    for i in range(1, N_QUEUE):
        ss[i + 1], xs[i + 1], ys[i + 1], txs[i + 1], tys[i + 1] = next_point(
            ss[i], xs[i], ys[i], L1, b, s_max
        )

    vs = [0.0] * (N_QUEUE + 1)
    vs[0] = V_HEAD
    for i in range(N_QUEUE):
        vs[i + 1] = (
            vs[i]
            * (txs[i] * (xs[i + 1] - xs[i]) + tys[i] * (ys[i + 1] - ys[i]))
            / (txs[i + 1] * (xs[i + 1] - xs[i]) + tys[i + 1] * (ys[i + 1] - ys[i]))
        )
    return xs, ys, vs


@jit
def min_distance(r, d):
    def point_to_segment_distance(p, a, b):
        # 计算点p到线段ab的最短距离，p在ab外侧时标记为inf
        px, py = p
        ax, ay = a
        bx, by = b

        abx = bx - ax
        aby = by - ay
        apx = px - ax
        apy = py - ay

        dot_product = apx * abx + apy * aby
        if dot_product <= 0:
            return float("inf")

        ab_len_sq = abx * abx + aby * aby
        t = dot_product / ab_len_sq
        if t >= 1:
            return float("inf")

        proj_x = ax + t * abx
        proj_y = ay + t * aby

        return sqrt((px - proj_x) ** 2 + (py - proj_y) ** 2)

    def vertices_to_polyline_distances(verts, polyline):
        # 计算顶点数组中每个点到折线的最短距离
        distances = []
        for v in verts:
            min_dist = float("inf")
            for i in range(len(polyline) - 1):
                a = polyline[i]
                b = polyline[i + 1]
                dist = point_to_segment_distance(v, a, b)
                if dist < min_dist:
                    min_dist = dist
            distances.append(min_dist)
        return np.array(distances)

    xs, ys, _ = queue(r, d)
    p0_x, p0_y = xs[0], ys[0]  # 龙头坐标
    p1_x, p1_y = xs[1], ys[1]  # 龙头后第一个把手
    p2_x, p2_y = xs[2], ys[2]  # 龙头后第二个把手
    vec10_x, vec10_y = p0_x - p1_x, p0_y - p1_y  # 龙头与龙头后第一个把手的连线方向
    vec21_x, vec21_y = p2_x - p1_x, p2_y - p1_y  # 第一个把手与第二个把手的连线方向

    vec_q_1 = sqrt(vec10_x**2 + vec10_y**2)
    vec_q_2 = sqrt(vec21_x**2 + vec21_y**2)
    n1_x, n1_y = vec10_x / vec_q_1, vec10_y / vec_q_1  # 单位化
    n2_x, n2_y = vec21_x / vec_q_2, vec21_y / vec_q_2

    vl0_x = p0_x + n1_x * H1 - n1_y * H2
    vl0_y = p0_y + n1_y * H1 + n1_x * H2  # 龙头左前方顶点位置
    vr0_x = p0_x + n1_x * H1 + n1_y * H2
    vr0_y = p0_y + n1_y * H1 - n1_x * H2  # 龙头右前方顶点位置

    vlf1_x = p1_x + n2_x * H1 - n2_y * H2
    vlf1_y = p1_y + n2_y * H1 + n2_x * H2  # 龙头后第一个把手左前顶点位置
    vrf1_x = p1_x + n2_x * H1 + n2_y * H2
    vrf1_y = p1_y + n2_y * H1 - n2_x * H2  # 龙头后第一个把手右前顶点位置

    vlr1_x = p1_x - n1_x * H1 - n1_y * H2
    vlr1_y = p1_y - n1_y * H1 + n1_x * H2  # 龙头后第一个把手左后顶点位置
    vrr1_x = p1_x - n1_x * H1 + n1_y * H2
    vrr1_y = p1_y - n1_y * H1 - n1_x * H2  # 龙头后第一个把手右后顶点位置

    verts = np.array(
        [
            [vl0_x, vl0_y],
            [vr0_x, vr0_y],
            [vlf1_x, vlf1_y],
            [vrf1_x, vrf1_y],
            [vlr1_x, vlr1_y],
            [vrr1_x, vrr1_y],
        ]
    )
    return min(vertices_to_polyline_distances(verts, np.array([xs[3:], ys[3:]]).T)) - H2


@jit
def d_r(r):
    dl, dr = 0.44, 0.46
    while dr - dl > 1e-12:
        dm = (dl + dr) / 2
        if min_distance(r, dm) * min_distance(r, dl) < 0:
            dr = dm
        else:
            dl = dm
    return (dl + dr) / 2


def find_max_valid_d():
    a, b = 4.5, 4.55
    gr = (sqrt(5) - 1) / 2
    c = b - (b - a) * gr
    d = a + (b - a) * gr
    while abs(b - a) > 1e-12:
        print(f"当前极径范围: [{a:.15f}, {b:.15f}]", end="\r")
        if d_r(c) > d_r(d):
            b = d
        else:
            a = c
        c = b - (b - a) * gr
        d = a + (b - a) * gr
    print()
    r = (a + b) / 2
    return r, d_r(r)


def plot_min_distance(rl, rr, dl, dr, nums):
    xs = np.linspace(rl, rr, nums)
    ys = np.linspace(dl, dr, nums)
    xss, yss = np.meshgrid(xs, ys)
    total = nums**2
    zs = np.zeros_like(xss)
    for i in range(nums):
        for j in range(nums):
            zs[i, j] = min_distance(xss[i, j], yss[i, j])
            print(f"绘制contour图: {i * nums + j + 1}/{total}", end="\r")
    print()
    fig = plt.figure(figsize=(10, 7))
    ax = fig.add_subplot(111)
    contourf = ax.contourf(
        xss, yss, zs,
        levels=20,
        cmap="viridis",
        alpha=0.7
    )
    plt.colorbar(contourf, ax=ax, label="最小距离")
    contour = ax.contour(
        xss, yss, zs,
        levels=20,
        colors="gray",
        linewidths=0.5
    )
    zero_contour = ax.contour(
        xss, yss, zs,
        levels=[0],
        colors="red",
        linewidths=2,
        linestyles="solid"
    )

    ax.set_title("最小距离与龙头极径、螺距的关系图")
    ax.set_xlabel("龙头极径")
    ax.set_ylabel("螺距")
    ax.grid(True, linestyle="--", alpha=0.5)
    ax.plot([4.5, 4.5], [dl, dr], color="red", linestyle="--")
    plt.title("最小距离与龙头极径、螺距的关系图")
    plt.savefig(f"p3-最小距离与龙头极径和螺距的关系图.svg")
    plt.clf()


plot_min_distance(4, 5, 0.4, 0.5, 80)
r_final, d_final = find_max_valid_d()
print(f"最小螺距: {d_final:.6f}，碰撞时龙头极径: {r_final:.6f}")