import numpy as np

N_M, N_FY, N_S = 3, 5, 3

ACC = 1e-2  # 计算精度
DT = ACC  # 时间积分步长
GOLDEN_RATIO = (np.sqrt(5) - 1) / 2  # 黄金分割比

G_GRAVITY = 9.8  # 重力加速度
G_GRAVITY_VECTOR = np.array((0, 0, -G_GRAVITY))

# 无人机初始位置
LOC_FY = np.array(
    (
        (17800, 0, 1800),
        (12000, 1400, 1400),
        (6000, -3000, 700),
        (11000, 2000, 1800),
        (13000, -2000, 1300),
    )
)
V_FY_LIMIT = np.array((70, 120))  # 无人机速度范围

V_SMOKE = 3  # 烟雾下降速度
V_SMOKE_VECTOR = np.array((0, 0, -V_SMOKE))
R_SMOKE = 10  # 烟雾半径
T_SMOKE = 20  # 烟雾存在时间

# 导弹初始位置与速度
LOC_M = np.array(
    (
        (20000, 0, 2000),
        (19000, 600, 2100),
        (18000, -600, 1900),
    )
)
LOC_M_MODULUS = np.zeros(N_M, dtype="float")
V_M = 300
V_M_VECTOR = np.zeros_like(LOC_M, dtype="float")
for _k in range(3):
    LOC_M_MODULUS[_k] = np.linalg.norm(LOC_M[_k])
    V_M_VECTOR[_k, :] = -LOC_M[_k] / LOC_M_MODULUS[_k] * V_M

LOC_REAL_TARGET = np.array((0, 200, 0))  # 真实目标位置
R_TARGET = 7  # 真实目标半径
H_TARGET = 10  # 真实目标高度

LOC_FAKE_TARGET = np.array((0, 0, 0))  # 诱饵目标位置

# 圆柱体取点个数
N_POINTS_ON_CIRS = 12
N_POINTS_ON_PART_CIR = N_POINTS_ON_CIRS // 4
N_POINTS_ON_CIR = N_POINTS_ON_PART_CIR * 2
N_POINTS_ON_CIRS = N_POINTS_ON_PART_CIR * 4
# 在圆柱的上下圆周的两侧取点
_D_ANGLE = np.pi / 18  # 取 10 度小角
_THETAS1 = np.linspace(
    1 / 2 * np.pi - _D_ANGLE, 1 / 2 * np.pi + _D_ANGLE, N_POINTS_ON_PART_CIR
)
_THETAS2 = np.linspace(
    3 / 2 * np.pi - _D_ANGLE, 3 / 2 * np.pi + _D_ANGLE, N_POINTS_ON_PART_CIR
)
_THETAS = np.hstack((_THETAS1, _THETAS2))
_THETAS.resize((N_POINTS_ON_CIR, 1))
_ONES = np.ones((N_POINTS_ON_CIR, 1))
# 下表面圆周
_POINTS_REAL_TARGET_DOWN = np.hstack(
    (R_TARGET * np.cos(_THETAS), R_TARGET * np.sin(_THETAS), 0 * _ONES)
)
# 上表面圆周
_POINTS_REAL_TARGET_UP = np.hstack(
    (R_TARGET * np.cos(_THETAS), R_TARGET * np.sin(_THETAS), H_TARGET * _ONES)
)
# 所有取点
POINTS_REAL_TARGET = np.concatenate((_POINTS_REAL_TARGET_DOWN, _POINTS_REAL_TARGET_UP))
POINTS_REAL_TARGET = POINTS_REAL_TARGET + LOC_REAL_TARGET

# 所有求解的最大起止时刻
T_START, T_END = 0, 67


# 预处理参数
R_SMOKE2 = R_SMOKE * R_SMOKE
V_M2 = V_M * V_M
