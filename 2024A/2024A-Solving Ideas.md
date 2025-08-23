# 问题1
## 数学模型
### 盘入螺线方程
设螺距为 $d$，采用极坐标形式，设极角为 $\theta$，极径为 $r$。由等距螺线的性质，极角每增加一圈，极径增加 $d$，二者呈线性关系。由初始条件，$\theta = 16 \cdot 2\pi$ 时，$r = 16d$，确定螺线方程为
$$
r = \frac{d}{2\pi} \theta
$$
### 直接递推关系（未使用）
设龙头前把手位置为 $P_0$，其后各把手依次为 $P_1, P_2, \cdots, P_N$，$N$ 为龙队节数。在极坐标系中，设 $P_i$ 的极角和极径分别为 $\theta_i$ 和 $r_i$。设龙队各节前后把手距离分别为 $l_1, l_2, \cdots, l_N$，则由余弦定理
$$
r_i^2 + r_{i+1}^2 - 2r_i r_{i+1} \cos(\theta_i - \theta_{i+1}) = l_{i+1}, \quad i = 0, 1, \cdots, N-1
$$

带入 $r = \frac{d}{2\pi} \theta$，整理得
$$
\theta_i^2 + \theta_{i+1}^2 - 2\theta_i \theta_{i+1} \cos(\theta_i - \theta_{i+1}) = 4\pi^2\frac{l_{i+1}^2}{d^2}, \quad i = 0, 1, \cdots, N-1
$$
在已知 $\theta_i$ 的情况下求解 $\theta_{i+1}$，该方程具有非常多零点。为获得正确零点，只需根据龙队运行方向，限定 $\theta_{i+1}$ 的范围为 $(\theta_i-\pi, \theta_i)$ 或 $(\theta_i, \theta_i+\pi)$。
由于问题4需要处理复合曲线，该递推关系式失效，实际计算中并未使用。需寻找更统一的求解方式。
### 螺线弧长公式
考虑螺线
$$
r = \frac{d}{2\pi} \theta
$$
记 $L(\theta)$ 为螺线上，极角在 $[0,\theta]$ 的弧长，则
$$
\begin{align}
L(\theta) 
&= \int_0^\theta \sqrt{r^2 + \left(\frac{\mathrm{d}r}{\mathrm{d}\theta}\right)^2} \mathrm{d}\theta \\
&= \frac{d}{2\pi} \int_0^\theta \sqrt{\theta^2 + 1}\mathrm{d}\theta \\
&= \frac{d}{4\pi} \left[\theta \sqrt{\theta^2 + 1} + \ln\left(\theta + \sqrt{\theta^2 + 1}\right)\right]
\end{align}
$$
$L$ 是 $\theta$ 的单调函数，存在反函数 $L^{-1}$。$L^{-1}$ 可通过牛顿迭代法快速求解。
### 旋转矩阵
为便于后续表述，设左乘后，将向量逆时针旋转 $\omega$ 的旋转矩阵为
$$
\mathbf{T}(\omega) =
\begin{bmatrix}
\cos(\omega) & -\sin(\omega) \\
\sin(\omega) & \cos(\omega) \\
\end{bmatrix}
$$
设 $x$，$y$ 轴方向的单位向量为 $\mathbf{e}_x$，$\mathbf{e}_y$。
旋转矩阵的导数为
$$
\frac{\mathrm{d}\mathbf{T}(\omega)}{\mathrm{d}\omega}
=
\begin{bmatrix}
-\sin(\omega) & -\cos(\omega) \\
\cos(\omega) & -\sin(\omega) \\
\end{bmatrix}
=
\mathbf{T}(\omega+\frac{\pi}{2})
$$
设 $\mathbf{a}$，$\mathbf{b}$ 为向量，旋转矩阵具有如下性质，
1. $\mathbf{T}(\omega)^{-1} = \mathbf{T}(-\omega)$；
2. $\mathbf{T}(\omega_1)\mathbf{T}(\omega_2) = \mathbf{T}(\omega_1+\omega_2)$；
3. $(\mathbf{T}(\omega_1)\mathbf{a})\cdot(\mathbf{T}(\omega_2)\mathbf{b})=(\mathbf{T}(\omega_1)\mathbf{a})^T(\mathbf{T}(\omega_2)\mathbf{b})=\mathbf{a}^T\mathbf{T}(\omega_1)^T\mathbf{T}(\omega_2)\mathbf{b}=\mathbf{a}^T\mathbf{T}(\omega_2-\omega_1)\mathbf{b}$

### 把手位置的自然坐标表示
记 $s_0$ 为龙头在螺线上的自然坐标，$s_1, s_2, \cdots, s_N$ 分别为龙头后各把手的自然坐标。以问题一中设定的龙头初始位置 $(r,\theta) = (16d,32\pi)$ 为自然坐标原点， 以龙队行进方向，即顺时针方向为自然坐标正方向。记 $\mathbf{r}_i$ 为第 $i$ 个把手的位置矢量。
设 $\theta_{max} = 32\pi$，为自然坐标起点在螺线上的极角。由自然坐标以曲线弧长为度量方式的性质，有
$$
s = L(\theta_{max}) - L(\theta)
$$
记 $L_{max} = L(\theta_{max}) = L(32\pi)$ 于是 $\theta$ 可表示为 $s$ 的函数
$$
\theta(s) = L^{-1}\left(L_{max} - s\right)
$$
$\theta(s)$ 的导数为
$$
\begin{align}
\frac{\mathrm{d}\theta(s)}{\mathrm{d}s}
&= \left(\frac{\mathrm{d}s}{\mathrm{d}\theta}\right)^{-1} \\
&= \left(\frac{\mathrm{d}}{\mathrm{d}\theta} \left[L(\theta_{max}) - L(\theta)\right]\right)^{-1} \\
&= \left(-\frac{\mathrm{d}L}{\mathrm{d}\theta} (\theta)\right)^{-1} \\
&= -\left(\frac{d}{2\pi} \sqrt{1 + \theta^2}\right)^{-1} \\
&= -\frac{2\pi}{d} \frac{1}{\sqrt{1 + \theta(s)^2}} \\
\end{align}
$$
记向量值函数 $\mathbf{r}(s)=\left[x,y\right]^T$ 表示：螺线上自然坐标为 $s$ 的点的位置矢量，则
$$
\mathbf{r}(s) = 
\begin{bmatrix}
x(s)\\
y(s)\\
\end{bmatrix}
=
\begin{bmatrix}
r(s)\cos(\theta(s))\\
r(s)\sin(\theta(s))\\
\end{bmatrix}
=
\frac{d}{2\pi} \theta(s)
\begin{bmatrix}
\cos(\theta(s)) & -\sin(\theta(s)) \\
\sin(\theta(s)) & \cos(\theta(s)) \\
\end{bmatrix}
\begin{bmatrix}
1 \\
0 \\
\end{bmatrix}
=
\frac{d}{2\pi} \theta(s) \mathbf{T}(\theta(s)) \mathbf{e}_x
$$
$\mathbf{r}(s)$ 的导数为
$$
\begin{align}
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
&= \frac{\mathrm{d}}{\mathrm{d}s}\left[\frac{d}{2\pi} \theta(s) \mathbf{T}(\theta(s)) \mathbf{e}_x\right] \\
&= \frac{d}{2\pi} \left[ \frac{\mathrm{d}\theta(s)}{\mathrm{d}s} \mathbf{T}(\theta(s)) + \theta(s)\frac{\mathrm{d}\mathbf{T}(\theta(s))}{\mathrm{d}s} \right]\mathbf{e}_x \\
&= \frac{d}{2\pi} \frac{\mathrm{d}\theta(s)}{\mathrm{d}s} \left[ \mathbf{T}(\theta(s)) + \theta(s)\mathbf{T}(\theta(s)+\frac{\pi}{2}) \right]\mathbf{e}_x \\
&= -\frac{1}{\sqrt{1 + \theta(s)^2}}\left[ \mathbf{T}(\theta(s)) + \theta(s)\mathbf{T}(\theta(s)+\frac{\pi}{2}) \right]\mathbf{e}_x
\end{align}
$$
展开 $x$，$y$ 分量计算，可得
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
=
-\frac{1}{\sqrt{1 + \theta(s)^2}}
\begin{bmatrix}
\cos(\theta(s)) - \theta(s) \sin(\theta(s)) \\
\sin(\theta(s)) + \theta(s) \cos(\theta(s)) \\
\end{bmatrix}
= -\frac{1}{\sqrt{1 + \theta(s)^2}} \mathbf{T}(\theta(s))
\begin{bmatrix}
1 \\
\theta(s) \\
\end{bmatrix}
$$
令 $\mathbf{b}(x) = \frac{1}{\sqrt{1 + x^2}} [1, x]^T$，则 $\mathbf{b}(x)$ 恰好为单位向量，且上式可简记为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s) = - \mathbf{T}(\theta(s)) \mathbf{b}(\theta(s))
$$
$\mathbf{b}(x)$ 也可进一步统一写为
$$
\mathbf{b}(x) = \mathbf{T}(\arctan(x))\mathbf{e}_x
$$
但这种表示方式会降低计算效率，实际并未使用。
至此，已经完全用自然坐标表示了把手的位置矢量。这种表示方法能极大提高计算效率，并可以直接适配问题4、5的情景。
### 基于自然坐标的递推关系
**方程推导**
在已知前一个把手的自然坐标 $s_i$ 的情况下，为了求解下一个把手的自然坐标 $s_{i+1}$ ，只需求解方程
$$
\|\mathbf{r}(s_{i+1}) - \mathbf{r}(s_i)\| = l_{i+1}
$$
其中，$l_i(1 \le i \le N)$ 为第 $i$ 条板凳前后把手之间的距离。 
该方程同样有众多零点。本问题中，自然坐标系的定义决定了必然有 $s_{i+1} < s_i$，因此只需求解小于 $s_i$ 的最大零点。设龙队各节前后把手距离分别为 $l_1, l_2, \cdots, l_N$，递推关系可表示为
$$
s_{i+1} = \max \set{s|\ \|\mathbf{r}(s) - \mathbf{r}(s_i)\| = l_{i+1}\ \text{且}\ s<s_i}
$$
**方程求解**
为求解方程，作如下辅助函数
$$
f(s,s_i,l_{i+1})
= \|\mathbf{r}(s) - \mathbf{r}(s_i)\|^2 - l_{i+1}^2
= [\mathbf{r}(s) - \mathbf{r}(s_i)]^2 - l_{i+1}^2
$$
原方程与 $f$ 具有相同的零点。
$f$ 对 $s$ 的导数为
$$
\frac{\mathrm{d}f}{\mathrm{d}s} = 2[\mathbf{r}(s) - \mathbf{r}(s_i)] \cdot \frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
$$
可使用牛顿迭代法求解，选取 $s_{i+1}$ 的迭代初值为 $s_i-\pi l_{i+1}$。
### 模型总结
令 $d$ 为螺距；
令 $v$ 为龙头速度；
令 $l_i$ 为各节龙身前后把手距离；
令 $L(\theta) = \frac{d}{4\pi} \left[\theta \sqrt{\theta^2 + 1} + \ln\left(\theta + \sqrt{\theta^2 + 1}\right)\right]$；
令 $\theta(s) = L^{-1}\left(L_{max} - s\right)$，其中 $L_{max} = L(32\pi)$；
位置矢量关于自然坐标的函数为
$$
\mathbf{r}(s)
= \frac{d}{2\pi} \theta(s) \mathbf{T}(\theta(s)) \mathbf{e}_x
$$
令 $\mathbf{b}(x) = \frac{1}{\sqrt{1 + x^2}} [1, x]^T$，位置矢量关于自然坐标的导数为
$$

\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s) = -\mathbf{T}(\theta(s)) \mathbf{b}(\theta(s))
$$
令
$$
f(s,s_i,l_{i+1})
= [\mathbf{r}(s) - \mathbf{r}(s_i)]^2 - l_{i+1}^2
$$
$f$ 对 $s$ 的导数为
$$
\frac{\mathrm{d}f}{\mathrm{d}s} = 2[\mathbf{r}(s) - \mathbf{r}(s_i)] \cdot \frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
$$
各把手自然坐标的递推式为
$$
\begin{cases}
s_0 = vt \\
s_{i+1} = \max \set{s|f(s,s_i,l_{i+1})=0\ \text{且}\ s<s_i}, \ i = 0, 1, \cdots, N
\end{cases}
$$
上述递推方程可通过牛顿迭代法求解，选取 $s_{i+1}$ 的迭代初值为 $s_i-\pi l_{i+1}$，即可完全保证得到需要的零点。
# 问题2
## 数学模型
### 板凳顶点位置
设板凳前后把手到前后沿的距离为 $h_1$，设把手到板凳左右边缘的距离为 $h_2$。
以龙头为第 $1$ 条板凳，前后把手位置矢量分别为 $\mathbf{r}(s_0)$，$\mathbf{r}(s_1)$。以此类推，第 $i(1 \le i \le N)$ 条板凳的前后把手的位置矢量分别为 $\mathbf{r}(s_{i-1})$，$\mathbf{r}(s_i)$。
第 $i$ 条板凳的方向为 $\mathbf{r}(s_i) - \mathbf{r}(s_{i+1})$，设 $\mathbf{n}_i$ 为指向第 $i$ 条板凳正前方的单位向量，则
$$
\mathbf{n}_i
= \frac{\mathbf{r}(s_i) - \mathbf{r}(s_{i+1})}{\left\|\mathbf{r}(s_i) - \mathbf{r}(s_{i+1})\right\|}
= \frac{1}{l_i}[\mathbf{r}(s_i) - \mathbf{r}(s_{i+1})]
$$
为描述板凳的四个顶点，还需求板凳的正右侧方向，令
$$
\mathbf{m}_i = \mathbf{T}(-\frac{\pi}{2})\mathbf{n}_i
$$
则 $\mathbf{m}_i$ 为 $\mathbf{n}_i$ 沿顺时针旋转 $\frac{\pi}{2}$ 后的单位向量，即指向第 $i$ 条板凳的正右方的单位向量。
第 $i$ 条板凳的四个顶点位置矢量分别为
$$
\begin{align}
\mathbf{p}_{i,\text{左前}} &= \mathbf{r}(s_{i-1}) + h_1\mathbf{n}_i - h_2\mathbf{m}_i \\
\mathbf{p}_{i,\text{右前}} &= \mathbf{r}(s_{i-1}) + h_1\mathbf{n}_i + h_2\mathbf{m}_i \\
\mathbf{p}_{i,\text{左后}} &= \mathbf{r}(s_i) - h_1\mathbf{n}_i - h_2\mathbf{m}_i \\
\mathbf{p}_{i,\text{右后}} &= \mathbf{r}(s_i) - h_1\mathbf{n}_i + h_2\mathbf{m}_i \\
\end{align}
$$
### 碰撞判断
由于各节龙身的长度相同，在前进过程中，各节龙身都会经历和第一节龙身相同的过程，无需重复判断龙身是否会发生碰撞。
对与长方体之间的碰撞，由于互不平行，碰撞形式必然是定点与边的碰撞。
又因为第一节龙身的后顶点和龙头的后顶点会经历相同的过程，故最终只需考虑龙头的前后顶点、第一节龙身的前顶点。由于右侧顶点靠内，盘入时不会发生碰撞，故只需考虑左侧的前三个顶点： $\mathbf{p}_{1,\text{左前}}$、$\mathbf{p}_{1,\text{左后}}$、$\mathbf{p}_{2,\text{左前}}$，下记为 $\mathbf{p}_1$、$\mathbf{p}_2$、$\mathbf{p}_3$。
经过以上分析，只需考虑三个顶点与第二节龙身及之后的板凳的位置关系。
对第 $i$ 条板凳，记板凳长向轴线分别与板凳的前、后边交点的位置矢量分别为 $\mathbf{a}_i$，$\mathbf{b}_i$。
则
$$
\begin{align}
\mathbf{a}_i &= \mathbf{r}(s_{i-1}) + h_1\mathbf{n}_i \\
\mathbf{a}_i &= \mathbf{r}(s_i) - h_1\mathbf{n}_i \\
\end{align}
$$
由于前两节板凳不会相互碰撞，只需考虑 $i \ge 3$ 的情况。
以 $c_{ij}$ 表示 $\mathbf{p}_j(j=1,2,3)$ 与第 $i(i \ge 3)$ 条板凳的距离。
当 $(\mathbf{p}_j - \mathbf{a}_i) \cdot (\mathbf{b}_i - \mathbf{a}_i) \le 0$ 或 $(\mathbf{p}_j - \mathbf{b}_i) \cdot (\mathbf{a}_i - \mathbf{b}_i) \le 0$ 时，点 $\mathbf{p}_j$ 到 $\mathbf{a}_i$，$\mathbf{b}_i$ 连线的垂足在 $\mathbf{a}_i$，$\mathbf{b}_i$ 的外部，此时必然不会发生碰撞，令 $c_{ij} = + \infty$。
当 $(\mathbf{p}_j - \mathbf{a}_i) \cdot (\mathbf{b}_i - \mathbf{a}_i) \ge 0$ 且 $(\mathbf{p}_j - \mathbf{b}_i) \cdot (\mathbf{a}_i - \mathbf{b}_i) \ge 0$ 时，点 $\mathbf{p}_j$ 到 $\mathbf{a}_i$，$\mathbf{b}_i$ 连线的垂足在 $\mathbf{a}_i$，$\mathbf{b}_i$ 的内部，距离 $c_{ij}$ 可表示为
$$
c_{ij} = \frac{\left\| (\mathbf{p}_j - \mathbf{a}_i) \times (\mathbf{b}_i - \mathbf{a}_i)  \right\|}{\left\| (\mathbf{b}_i - \mathbf{a}_i)  \right\|} - h_2
$$
发生碰撞时，存在 $c_{ij} \le 0$，令
$$
c_{min} = \min\set{c_{ij}|1 \le i \le N, 1 \le j \le 3}
$$
则发生的时刻对应了 $c_{min}=0$ 的时刻。
由问题1中的递推关系式 $c_{min}$ 以龙头坐标 $s_0$ 为唯一变量。又因为 $s_0 = vt$，$c_{min}$ 也是 $t$ 的函数，可写作  $c_{min}(t)$。记第一次碰撞的时刻为 $t_0$，则 $t_0$ 为方程 $c_{min}(t)=0$ 最小的零点，于是
$$
t_0 = \min\set{t|c_{min}(t)=0}
$$
# 问题3
## 数学模型
### 约束条件
设掉头区域半径为 $R$。
对于任意螺距 $d$，龙头的极径 $c_{min}$ 可视作时刻 $t$ 与螺距 $d$ 的函数，由于每个 $t$ 对应了唯一的龙头极径 $r_0$，故 $c_{min}$ 可也表示为 $r_0$ 与 $d$ 的函数，记为 $c_{min}(r_0, d)$。
$r_0$ 与 $t$ 的关系如下
$$
r_0 = \frac{d}{2\pi} \theta_0 = \frac{d}{2\pi} L^{-1}\left(L_{max} - vt\right)
$$
$$
t = \frac{L_{max} - \frac{2\pi r_0}{d}}{v}
$$
为保证在龙头到达掉头区域前，队伍不发生碰撞，$d$ 需满足：当$r_0 \ge R$ 时，$c_{min}(r_0, d)$ 恒大于 $0$。
### 求解过程
绘制函数 $c_{min}(r_0, d)$，的着色图，如下所示
![[p3-最小距离与龙头极径和螺距的关系图.svg]]
在实际运动时，螺距为固定值，极径将由大变小。上图中，表示队伍状态的点将从右至左沿水平线运动。图中红实线表示 $c_{min}(r_0, d)=0$ 的等高线，红色虚线表示 $r_0 = R$，为使碰撞不发生，表示队伍状态的点在从右侧水平运动到红色虚线的过程中，必须保持在红色实线上方。
上图中的红色实线对应了一个 $d$ 关于 $r_0$ 的单值函数，记为 $D(r_0)$，注意到在红色虚线右侧附近 $D(r_0)$ 为上凸函数，容易求得极大值 $D_{max}$，此极大值即为最小允许的螺距 $d_{min}$。
# 问题4
## 数学模型
### 弧线长与半径比例的关系
记掉头开始和结束的点分别为 $A$，$B$，圆弧交界点为 $M$，记两圆的圆心分别为 $O_1$，$O_2$。作 $OH$ 垂直点 $A$ 处的切线，交切线于点 $H$。
如下图所示
![[p4-示意图.svg]]
为简单证明弧长不变，记 $R_1 = O_1M$，$R_2 = O_2M$，有
$$
\begin{align}
BH &= (R_1 + R_2) + (R_1 + R_2)\cos(\angle HO_2O_1)\\
AH &= (R_1 + R_2)\sin(\angle HO_2O_1)
\end{align}
$$
解得
$$
\begin{align}
R_1 + R_2 &= \frac{AH^2+BH^2}{2BH} \\
\angle HO_2O_1 &= 2 \arctan{\frac{AH}{BH}}
\end{align}
$$
两段圆弧长度的和为
$$
C = 2(R_1 + R_2)(\pi - \angle HO_2O_1) = \frac{AH^2+BH^2}{2BH}\left(\pi -2 \arctan{\frac{AH}{BH}}\right)
$$
为定值。
### 曲线的拼接
队伍的轨迹由四段不同的曲线拼接而成，现需要重新定义运动轨迹的自然坐标。以点 $A$ 为运动曲线自然坐标的坐标原点，以队伍行进方向为自然坐标正方向。整条运动曲线由四部分拼接而成，以下分别讨论。
**第一段曲线**
第一段曲线上，$s \le 0$，把手在螺线上运动。
掉头区域半径为 $R$，记螺线和掉头区域圆周相交的点为 $A$，记由原点 $O$ 到 $A$ 的螺线的弧长为 $L_{min}$，记点 $A$ 的极角和极径为 $\theta_{min}$、$r_{min}$，则
$$
\begin{align}
r_{min} &= R \\
\theta_{min} &= \frac{2\pi}{d} r_{min} = \frac{2\pi}{d} R \\
L_{min} &= L(\theta_{min})
\end{align}
$$
由于 $A$ 为自然坐标原点，$A$ 的位置矢量为 $\mathbf{r}(0)$。
记把手的角和极径分别为 $\theta$，$r$，则在第一段曲线上，$\theta$ 可视为 $s$ 的函数，记作 $\theta(s)$。
由
$$
s = -(L(\theta)-L_{min}) = L_{min} - L(\theta)
$$
有
$$
\theta(s) = L^{-1}(L_{min} - s)
$$
使用重新定义的 $\theta(s)$，与问题1中相同，把手的位置矢量表示为
$$
\mathbf{r}(s)
= \frac{d}{2\pi} \theta(s) \mathbf{T}(\theta(s)) \mathbf{e}_x
,\quad s \le 0
$$
单位切向量为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
= - \mathbf{T}(\theta(s)) \mathbf{b}(\theta(s))
,\quad s \le 0
$$
其中，$\mathbf{b}(x) = \frac{1}{\sqrt{1+x^2}}[1, x]^T$。
**第二、三段曲线**
第二、三段曲线均为圆弧，设两条圆弧所在圆的圆心分别为 $O_1$，$O_2$，半径分别为 $R_1$，$R_2$，且令
$$
R_1:R_2 = \lambda_1:\lambda_2,\quad \lambda_1 + \lambda_2 = 1
$$
在问题4的情景中，$\lambda_1 = \frac{2}{3}$，$\lambda_2 = \frac{1}{3}$。
第一段圆弧的起点为 $A$，$A$ 点的自然坐标为 $0$，位置矢量为
$$
\mathbf{r}(0)
= \frac{d}{2\pi} \theta(0) \mathbf{T}(\theta(0)) \mathbf{e}_x
= \frac{d}{2\pi} \theta_{min} \mathbf{T}(\theta_{min}) \mathbf{e}_x
= R \mathbf{T}(\theta_{min}) \mathbf{e}_x
$$
令 $\mathbf{n} = \mathbf{T}(\theta_{min}) \mathbf{e}_x$，$\mathbf{n}$ 为单位向量，且
$$
\overrightarrow{OA} = \mathbf{r}(0) = R\mathbf{n}
$$
$A$ 点的切向单位矢量为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(0)
= - \mathbf{T}(\theta(0)) \mathbf{b}(0)
= - \mathbf{T}(\theta_{min}) \mathbf{b}(\theta_{min})
$$
第二段圆弧的终点为 $B$，$B$ 与 $A$ 关于原点对称，位置矢量相反、切向矢量相同。
$\overrightarrow{AO_1}$ 指向 $A$ 点切向矢量的正右侧，于是 $\overrightarrow{AO_1}$ 方向的单位向量为
$$
\mathbf{m}
= \mathbf{T}\left(-\frac{\pi}{2}\right) \frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(0)
= -\mathbf{T}\left(\theta_{min}-\frac{\pi}{2}\right) \mathbf{b}(\theta_{min})
= \mathbf{T}\left(\theta_{min}+\frac{\pi}{2}\right) \mathbf{b}(\theta_{min})
$$
两圆的半径分别为 $R_1$，$R_2$，于是
$$
\begin{align}
\overrightarrow{AO_1} &= R_1 \mathbf{m} \\
\overrightarrow{BO_2} &= -R_2 \mathbf{m} \\
\end{align}
$$
由此可得
$$
\begin{align}
\left\| \overrightarrow{O_2O_1} \right\|^2
&= \left\| \overrightarrow{AO_1} - \overrightarrow{BO_2} - \overrightarrow{AB}\right\|^2 \\
&= (R_1 + R_2)^2 \mathbf{m}^2 - 2(R_1 + R_2) \mathbf{m} \cdot \overrightarrow{AB} + \overrightarrow{AB} \cdot \overrightarrow{AB} \\
&= (R_1 + R_2)^2 - 2(R_1 + R_2) \mathbf{m} \cdot \overrightarrow{AB} + \overrightarrow{AB} \cdot \overrightarrow{AB}
\end{align}
$$
因为两圆弧相切，$\|O_1O_2\| = R_1 + R_2$，于是
$$
(R_1 + R_2)^2 = (R_1 + R_2)^2 - 2(R_1 + R_2) \mathbf{m} \cdot \overrightarrow{AB} + \overrightarrow{AB} \cdot \overrightarrow{AB}
$$
解得
$$
R_1 + R_2 = \frac{\overrightarrow{AB} \cdot \overrightarrow{AB}}{2 \mathbf{m} \cdot \overrightarrow{AB}}
$$
由 $A$，$B$ 关于原点对称 $\overrightarrow{AB} = 2 \overrightarrow{AO} = -2R\mathbf{n}$，带入化简可得
$$
\begin{align}
R_1 + R_2
&= \frac{(2R\mathbf{n})^2}{-2\mathbf{m} \cdot 2R\mathbf{n}} \\
&= -\frac{R}{[\mathbf{T}\left(\theta_{min}+\frac{\pi}{2}\right) \mathbf{b}(\theta_{min})] \cdot [\mathbf{T}(\theta_{min}) \mathbf{e}_x]} \\
&= -\frac{R}{\mathbf{b}(\theta_{min})^T \mathbf{T}\left(-\theta_{min}-\frac{\pi}{2}\right) \mathbf{T}(\theta_{min}) \mathbf{e}_x} \\
&= \frac{R}{\mathbf{b}(\theta_{min})^T \mathbf{e}_y} \\
&= \sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}
\end{align}
$$
由 $R_1 : R_2 = \lambda_1:\lambda_2$，且满足 $\lambda_1+\lambda_2=1$，有
$$
\begin{align}
R_1 &= \frac{\lambda_1}{\mathbf{b}(\theta_{min}) \mathbf{e}_y} R
= \lambda_1\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2} \\
R_2 &= \frac{\lambda_2}{\mathbf{b}(\theta_{min}) \mathbf{e}_y} R 
= \lambda_2\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}\\
\end{align}
$$
综合以上等式，可得两圆圆心和交点的位置矢量分别为
$$
\begin{align}
\overrightarrow{OO_1} &= \overrightarrow{OA} + \overrightarrow{AO_1} = R \mathbf{n} + R_1 \mathbf{m} \\
\overrightarrow{OO_2} &= \overrightarrow{OB} + \overrightarrow{BO_2} = -R \mathbf{n} - R_2 \mathbf{m} \\
\overrightarrow{OM} &= \lambda_2 \overrightarrow{OO_1} + \lambda_1 \overrightarrow{OO_2} = (\lambda_2-\lambda_1) R\mathbf{n} \\
\end{align}
$$
其中
$$
\begin{align}
\mathbf{n} &= \mathbf{T}(\theta_{min}) \mathbf{e}_x \\
\mathbf{m} &= \mathbf{T}\left(\theta_{min}+\frac{\pi}{2}\right) \mathbf{b}(\theta_{min}) \\
\end{align}
$$
由几何平行关系，两端圆弧圆心角相等，均为
$$
\beta = \angle MO_1A = 2\arctan\left( \frac{\overrightarrow{O_1M} \times \overrightarrow{O_1A}}{\overrightarrow{O_1M} \cdot \overrightarrow{O_1A}} \right)
$$
化简得到
$$
\beta = 2\arctan{\theta_{min}} = 2\arctan\left(\frac{2\pi R}{d}\right)
$$
两段圆弧长分别为
$$
C_1 = \beta R_1 = 2\lambda_1\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}\arctan\left(\frac{2\pi R}{d}\right)
$$
$$
C_2 = \beta R_2 = 2\lambda_2\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}\arctan\left(\frac{2\pi R}{d}\right)
$$
圆弧长度之和为
$$
C = C_1 + C_2 = 2\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}\arctan\left(\frac{2\pi R}{d}\right)
$$
在第二段曲线上，$0 < s \le C_1$，把手在圆弧上顺时针运动，走过的圆心角为 $\frac{s}{R_1}$，位置矢量 $\mathbf{r}$ 可表示为
$$
\begin{align}
\mathbf{r}
&= \overrightarrow{OO_1} + \mathbf{T}\left(-\frac{s}{R_1}\right) \overrightarrow{O_1A} \\
&= R \mathbf{n} + R_1 \mathbf{m} - \mathbf{T}\left(-\frac{s}{R_1}\right) R_1 \mathbf{m} \\
\end{align}
$$
单位切向量为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
= \mathbf{T}\left(\frac{\pi}{2}-\frac{s}{R_1}\right) \mathbf{m} \\
$$
在第三段曲线上，$C_1 < s < C_1+C_2$，把手在圆弧上逆时针运动，距离 $B$ 点的圆心角为 $\frac{C_1 + C_2 - s}{R_1}$，同理，位置矢量 $\mathbf{r}$ 可表示为
$$
\begin{align}
\mathbf{r}(s)
&= \overrightarrow{OO_2} + \mathbf{T}\left(-\frac{C_1+C_2-s}{R_2}\right) \overrightarrow{O_2B} \\
&= -R \mathbf{n} - R_2 \mathbf{m} + \mathbf{T}\left(-\frac{C_1+C_2-s}{R_2}\right) R_2 \mathbf{m}
\end{align}
$$
单位切向量为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
= -\mathbf{T}\left(\frac{\pi}{2}-\frac{C_1+C_2-s}{R_2}\right) \mathbf{m} \\
$$
**第四段曲线**
在第四段曲线上，$s \ge C_1+C_2$，把手由内至外绕盘出螺线运动，令 $s'=C_1+C_2-s$，与第一段曲线同理，可得
$$
\mathbf{r}(s)
= -\frac{d}{2\pi} \theta(C_1+C_2-s) \mathbf{T}(\theta(C_1+C_2-s)) \mathbf{e}_x
,\quad s \ge C_1+C_2
$$
单位切向量为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
= \mathbf{T}(\theta(C_1+C_2-s)) \mathbf{b}(\theta(C_1+C_2-s))
,\quad s \ge C_1+C_2
$$
至此，完成了所有曲线信息的计算：
- 第一段曲线为盘入螺线，对应的自然坐标范围为 $(-\infty,0]$；
- 第二段曲线为以 $O_1$ 为圆心的圆弧 $AM$，对应的自然坐标范围为 $(0,C_1)$；
- 第三段圆弧为以 $O_2$ 为圆心的圆弧 $MB$，对应的自然坐标范围为 $[C_1,C_2)$；
- 第四段圆弧为盘出螺线，对应的自然坐标范围为 $[C_2,+\infty)$。

### 递推关系
同问题1，有递推关系
$$
\begin{cases}
s_0 = vt \\
s_{i+1} = \max \set{s|f(s,s_i,l_{i+1})=0\ \text{且}\ s<s_i}, \ i = 0, 1, \cdots, N
\end{cases}
$$
其中
$$
f(s,s_i,l_{i+1})
= [\mathbf{r}(s) - \mathbf{r}(s_i)]^2 - l_{i+1}^2
$$
由于曲线上各点的位置和单位切向量已知，且位置和单位切向量均连续，故同样可通过牛顿迭代法求解。
### 模型总结
令 $d$ 为螺距；令 $R$ 为掉头半径；令 $v$ 为龙头速度；
令 $l_i$ 为各节龙身前后把手距离；
令 $\lambda_1$，$\lambda_2$ 为两段圆弧半径的比例，满足 $\lambda_1 + \lambda_2 = 1$；
令 $R_1$，$R_2$ 分别为掉头区域两端圆弧的半径，则
$$
R_1 = \lambda_1\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2} ,\quad
R_2 = \lambda_2\sqrt{\left(\frac{d}{2\pi}\right)^2+R^2}
$$
令 $\beta$ 为两端圆弧的圆心角，则
$$
\beta = 2\arctan\left(\frac{2\pi R}{d}\right)
$$
令 $C_1$，$C_2$ 分别为两端圆弧的长度，则
$$
C_1 = \beta R_1 ,\quad
C_2 = \beta R_2
$$
令 $L(\theta) = \frac{d}{4\pi} \left[\theta \sqrt{\theta^2 + 1} + \ln\left(\theta + \sqrt{\theta^2 + 1}\right)\right]$；
令 $\theta_{min} = \frac{2\pi R}{d}$，$L_{min} = L(\theta_{min})$；
令 $\theta(s) = L^{-1}(L_{min} - s)$；
令 $\mathbf{T}(\omega)$ 为旋转 $\omega$ 的矩阵；
令 $\mathbf{b}(x) = \frac{1}{\sqrt{1 + x^2}} [1, x]^T$；
令 $\mathbf{n}$，$\mathbf{m}$ 分别为
$$
\mathbf{n} = \mathbf{T}(\theta_{min}) \mathbf{e}_x, \quad
\mathbf{m} = \mathbf{T}\left(\theta_{min}+\frac{\pi}{2}\right) \mathbf{b}(\theta_{min})
$$
则 $\mathbf{n}$，$\mathbf{m}$ 均为单位向量，且均为常量；
以开始进入掉头区域的位置为自然坐标原点，以龙队行进方向为自然坐标正方向，建立自然坐标系。
位置矢量关于自然坐标的函数为
$$
\mathbf{r}(s)
= 
\begin{cases}
\frac{d}{2\pi} \theta(s) \mathbf{T}(\theta(s)) \mathbf{e}_x &,s \le 0 \\
R \mathbf{n} + R_1 \mathbf{m} - R_1 \mathbf{T}\left(-\frac{s}{R_1}\right) \mathbf{m} &,0 \le s \le C_1 \\
-R \mathbf{n} - R_2 \mathbf{m} + R_2 \mathbf{T}\left(-\frac{C_1+C_2-s}{R_2}\right) \mathbf{m} &,C_1 \le s \le C_1 + C_2 \\
-\frac{d}{2\pi} \theta(C_1+C_2-s) \mathbf{T}(\theta(C_1+C_2-s)) \mathbf{e}_x&, C_1+C_2 \le s \\
\end{cases}
$$
位置矢量关于自然坐标的函数的导数为
$$
\frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
= 
\begin{cases}
- \mathbf{T}(\theta(s)) \mathbf{b}(\theta(s))
&,s \le 0 \\
\mathbf{T}\left(\frac{\pi}{2}-\frac{s}{R_1}\right) \mathbf{m}
&,0 \le s \le C_1 \\
-\mathbf{T}\left(\frac{\pi}{2}-\frac{C_1+C_2-s}{R_2}\right) \mathbf{m}
&,C_1 \le s \le C_1+C_2 \\
\mathbf{T}(\theta(C_1+C_2-s)) \mathbf{b}(\theta(C_1+C_2-s))
&, C_1+C_2 \le s \\
\end{cases}
$$
令
$$
f(s,s_i,l_{i+1})
= [\mathbf{r}(s) - \mathbf{r}(s_i)]^2 - l_{i+1}^2
$$
$f$ 对 $s$ 的导数为
$$
\frac{\mathrm{d}f}{\mathrm{d}s} = 2[\mathbf{r}(s) - \mathbf{r}(s_i)] \cdot \frac{\mathrm{d}\mathbf{r}}{\mathrm{d}s}(s)
$$
各把手自然坐标的递推式为
$$
\begin{cases}
s_0 = vt \\
s_{i+1} = \max \set{s|f(s,s_i,l_{i+1})=0\ \text{且}\ s<s_i}, \ i = 0, 1, \cdots, N
\end{cases}
$$
使用牛顿迭代法完成求解。计算 $s_{i+1}$ 时，选取的迭代初值为 $s_i-\pi l_{i+1}$。
# 问题5
## 数学模型
### 速度递推关系
问题4中的各个把手的自然坐标 $s_i(0 \le i \le N)$ 均可视为时间 $t$ 的函数。问题4中，已经求解了曲线所有部分的位置和位置关于自然坐标的导数，因此，可以设
$$
\begin{align}
&x_i = x(s_i), \  y_i = y(s_i) \\
&x_i' = \frac{\mathrm{d}x}{\mathrm{d}s}(s_i), \  y_i' = \frac{\mathrm{d}y}{\mathrm{d}s}(s_i) \\
&x_{i+1} = x(s_{i+1}), \  y_{i+1} = y(s_{i+1}) \\
&x_{i+1}' = \frac{\mathrm{d}x}{\mathrm{d}s}(s_{i+1}), \  y_{i+1}' = \frac{\mathrm{d}y}{\mathrm{d}s}(s_{i+1}) \\
\end{align}
$$
由等式
$$
(x(s_{i+1})-x(s_i))^2 + (y(s_{i+1})-y(s_i))^2 - l_{i+1}^2=0
$$
两端同时微分可得
$$
(x_{i+1}-x_i)(x_{i+1}'\mathrm{d}s_{i+1}-x_{i}'\mathrm{d}s_{i})
+ (y_{i+1}-y_i)(y_{i+1}'\mathrm{d}s_{i+1}-y_{i}'\mathrm{d}s_{i})
= 0
$$
整理得
$$
\frac{\mathrm{d}s_{i+1}}{\mathrm{d}s_i}
= -\frac{(x_{i+1}-x_i)x_{i}' + (y_{i+1}-y_i)y_{i}'}{(x_{i+1}-x_i)x_{i+1}' + (y_{i+1}-y_i)y_{i+1}'}
$$
由 $v_i = \frac{\mathrm{d}s_i}{\mathrm{d}t}$，$v_{i+1} = \frac{\mathrm{d}v_{i+1}}{\mathrm{d}t}$，有
$$
\frac{v_{i+1}}{v_i}
= \frac{\frac{\mathrm{d}s_{i+1}}{\mathrm{d}t}}{\frac{\mathrm{d}s_i}{\mathrm{d}t}}
= \frac{\mathrm{d}s_{i+1}}{\mathrm{d}s_i}
$$
设龙头速度为 $v$，于是有递推关系式
$$
\frac{v_{i+1}}{v_i}
= -\frac{(x_{i+1}-x_i)x_{i}' + (y_{i+1}-y_i)y_{i}'}{(x_{i+1}-x_i)x_{i+1}' + (y_{i+1}-y_i)y_{i+1}'}
, \ i = 0, 1, \cdots, N
$$
### 速度最大值描述
由递推关系式，速度的比值仅与各点的位置有关，而与时间无关，故队列中各把手的速度 $v_i(1 \le i \le N)$ 可视为龙头坐标 $s_0$ 与龙头速度 $v_0$ 的函数。记队列中最大的速度为 $v_{max}$，$v_{max}$ 也为  $s_0$ 和 $v_0$ 的函数，且有
$$
v_{max}(s_0,v_0) = \max\set{v_i(s_0,v_0)|1 \le i \le N}
$$
当 $s_0$ 固定时，$v_i(1 \le i \le N)$ 与 $v_0$ 的比值恒为定值，于是 $v_{max}$ 与 $v_0$ 的比值也恒为定值
$$
\frac{v_{max}(s_0,v_0)}{v_0} = \frac{v_{max}(s_0,1)}{1}
$$
$$
v_0 = \frac{v_{max}(s_0,v_0)}{v_{max}(s_0,1)}
$$
题目给出了 $v_{max}(s_0,v_0)$ 被允许的最大值，记为 $v_{quemax}$。为求 $v_0$ 被允许的最大值，只需计算 $v_{max}(s_0,1)$ 最大值，从而保证
$$
v_0 \max\set{v_{max}(s_0,1)|s_0 \in \mathbb{R}} \le v_{quemax}
$$
则 $v_0$ 的最大值
$$
v_{0max} = \frac{v_{quemax}}{\max\set{v_{max}(s_0,1)|s_0 \in \mathbb{R}}}
$$
### 龙头速度为 1 时，队列最大速度的计算
$v_{max}(s_0,1)$ 是仅与 $s_0$ 有关的函数，初步绘制其图像，如下图所示
![[p5-队列最大速度与龙头位置关系图.svg]]
在 $s_0 \ge 16$ 的范围内，图像呈近似的周期变化，且极值点逐步略微增大，但明显小于 $s_0 < 16$ 范围内的极值点。 
注意到该函数的极大值点在区间 $[14,15]$ 内，且函数在该区域内是上凸的，容易求得 $v_{max}(s_0,1)$ 的极大值。
由公式
$$
v_{0max} = \frac{v_{quemax}}{\max\set{v_{max}(s_0,1)|s_0 \in \mathbb{R}}}
$$
可计算被允许的龙头最大速度。