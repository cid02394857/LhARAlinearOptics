# LhARA 第一阶段束线线性光学模拟：代码解读与科研使用指南

> 适用于本仓库 `claude/busy-euler-mjur5i` 分支（2026-10）。
> 理论推导以 `00-Documentation/01-WriteUp/LhARA-TN-2024-07.pdf`（下称 **TN**，LaTeX 源码在同一目录）为准。
> 本指南按“物理 → 方法 → 代码 → 参数 → 运行 → 科研”的顺序讲解，并标出对应的 TN 章节，方便对照阅读。
> 文中所有命令都在 2026-10-23 实际运行验证过。

## 目录

1. [这套代码在模拟什么](#1-这套代码在模拟什么)
2. [方法论速览（TN 第 2–5 章）](#2-方法论速览tn-第-25-章)
3. [代码结构与数据流](#3-代码结构与数据流)
4. [参数文件详解（11-Parameters）](#4-参数文件详解11-parameters)
5. [运行：从生成粒子到出图](#5-运行从生成粒子到出图)
6. [如何用这套代码做有效的科学研究](#6-如何用这套代码做有效的科学研究)
7. [常见问题](#7-常见问题)
8. [附录：文件速查与 TN 章节对照](#8-附录文件速查与-tn-章节对照)

---

## 1. 这套代码在模拟什么

### 1.1 LhARA 与第一阶段

LhARA（Laser-hybrid Accelerator for Radiobiological Applications，用于放射生物学研究的激光混合加速器）的思路是：

1. 用高功率超短脉冲激光打击薄靶，通过**靶背法向鞘层加速（TNSA）**产生质子或离子。这种束流强度高、脉冲短，但发散角大、能谱很宽。
2. 用 **Gabor 透镜**（用磁场和电场约束电子等离子体、靠电子的空间电荷聚焦离子的透镜）俘获并聚焦这些粒子。
3. 经过能量选择和匹配，把束流送到实验终端站，做细胞照射等放射生物学实验。

- **第一阶段（Stage 1）**：激光源产生的质子（参数文件里的参考能量是 15 MeV）直接输运到**低能体外（in-vitro）终端站**。本代码模拟的就是这一段。
- **第二阶段**：用 FFA 环把束流加速到更高能量，再送往高能终端站。**不在本代码范围内。**

### 1.2 第一阶段束线的组成

以 `11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv` 为例（Section 是参数文件里的分段名）：

| Section | 主要元件 | 作用 |
|---|---|---|
| Source | 参数化 TNSA 源 | 产生粒子 |
| Interface | 漂移段，加两个圆形光阑（半径 2.0 mm、2.87 mm） | 靶室“喷嘴”，连接靶区和束线 |
| Capture | 2 个 Gabor 透镜（长 0.857 m） | 俘获大发散角的束流，把它变成近平行束 |
| Energy selection | 2 个 RF 腔（5 MV/m、200 MHz）、Gabor 透镜、准直器（半径 2.5 mm） | 利用色差（不同能量的粒子焦点位置不同），在准直器处选出所需能量 |
| Matching | 4 个 Gabor 透镜 | 把束流匹配到偏转段的入口 |
| Arc | 坐标切换（RPLCswitch）、2 块 45° 扇形二极铁、6 个四极铁 | 把束流偏转 90°，送往低能终端站 |
| 末端 | 2 m 漂移段 | 到达终端站 |

其他相关配置：

- `LhARABeamLine-Params-*-Solenoid.csv`：用螺线管代替 Gabor 透镜，方便和 Gabor 方案对比。
- `LhARABeamLine-Params-Gauss-*.csv`：源改成简单的高斯分布，便于检查束线本身的光学性能。
- `LhARA-Stage1-BeamDelivery.csv`：**从文件读入粒子**（例如 BDSIM 的输出），只模拟最后的“散焦 + 聚焦”四极铁对和到终端站的漂移。原作者 2026-09 的笔记显示，当前工作重点正是“低能终端站束流输运（beam delivery）的优化”，用的就是这一配置。

### 1.3 代码具体在算什么

1. **建束线**：按参数文件从上到下依次创建元件，让一个“参考粒子”（理想的 15 MeV 质子，沿设计轨道运动）先走一遍，确定每个元件入口和出口的位置、方向和动量。
2. **生成粒子**：从源抽样 N 个粒子（代码里称为 N 个“事件”）。能量取自 TNSA 能谱，方向和位置按模型分布抽样。
3. **追踪**：每个粒子的 6 维迹空间矢量逐个元件乘以 6×6 传输矩阵。四极铁、螺线管、二极铁和 Gabor 透镜的矩阵按**该粒子自己的动量**重新计算，所以能量不同的粒子受到的聚焦也不同（色差）。这正是能量选择的物理基础。
4. **丢失判断**：粒子撞上光阑、超出真空室半径，或者角度大到超出线性近似的适用范围，就判为丢失，不再追踪。
5. **记录与分析**：在每个元件出口记录全部存活粒子的迹空间坐标，写入二进制数据文件。分析程序再从数据文件计算束斑尺寸、发射度、Twiss 参数、传输效率等。

### 1.4 近似与局限（写论文时必须交代）

- **线性（一阶）光学、近轴近似**：不含高阶像差和边缘场。八极铁在线性近似下等价于一段漂移。
- **没有空间电荷效应**：激光源附近束流密度很高，这一段的结果需要谨慎对待。
- **Gabor 透镜按理想线性透镜处理**（电子密度均匀）。
- **光阑是“硬切割”**：不模拟粒子在材料里的散射和能损。
- **角度判据**：(x′² + y′²)/(2D²) > 1 的粒子会被丢弃，其中 D² = 1 + 2δ/β₀ + δ²（见 `BeamLineElement.ExpansionParameterFail`）。

所以，这套代码最适合做**快速设计、参数扫描、优化和交叉检验**。关键结论最好再用 BDSIM（基于 Geant4 的全追踪程序）确认。代码提供了 BDSIM 格式的导出（`2BDSIM.py`）和读入（`-z` 选项）。

---

## 2. 方法论速览（TN 第 2–5 章）

### 2.1 坐标系（TN §2）

- **实验室系**：原点在激光焦点（靶上）；z 轴水平，沿俘获轴指向下游；y 轴竖直向上；x 轴与 y、z 构成右手系。
- **参考粒子局部坐标系（RPLC）**：原点跟随参考粒子；z_r 沿参考粒子的运动方向，x_r 指向轨迹曲率中心的外侧。s 是参考粒子走过的路程。
- 元件依次放在参考粒子的轨迹上，需要时做坐标变换。`RPLCswitch` 元件用来在偏转前切换局部坐标轴的取向。

### 2.2 迹空间（TN §3）

每个粒子用 6 维迹空间矢量描述：

$$\phi = (x,\ x',\ y,\ y',\ z,\ \delta)$$

- x′ = p_x/p₀，y′ = p_y/p₀，即相对参考轨道的角度；
- z = Δs/β₀，即和参考粒子之间的纵向距离差；
- δ = ΔE/(c p₀)，即能量偏差。

这里 p₀、E₀、β₀ 是参考粒子在该位置的动量、能量和速度。

### 2.3 传输矩阵（TN §4）

粒子穿过一个元件：φ_出 = M φ_入，M 是 6×6 矩阵。代码里的元件类都在 `01-Code/BeamLineElement.py`：

| 类 | 参数文件里的写法（Element） | 关键参数 | 物理作用 |
|---|---|---|---|
| `Drift` | `Drift` | `Length` (m) | 自由漂移 |
| `FocusQuadrupole` / `DefocusQuadrupole` | `Fquad` / `Dquad` | `Length` (m)、`kq` (m⁻²) | 一个横向平面聚焦、另一个散焦。k_q = g/(Bρ)，g 为磁场梯度 (T/m)，Bρ 为磁刚度 |
| `Solenoid` | `Solenoid` | Type 为 `Length, strength` 时给 `Length`、`Strength` (rad/m)；Type 为 `Length, layers and turns` 时给 `Length`、`Layers`、`Turns`、`Current` (A) | 两个横向平面都聚焦，同时引起 x–y 耦合（转动）。`Strength` = B₀/(Bρ)；线圈写法用 B₀ = μ₀NI/L |
| `GaborLens` | `Gabor lens` | `Length` (m)、`Strength` (rad/m) | 两个平面都聚焦，没有转动。强度由电子密度决定：k_G = (e/2ε₀)·(m_p γ/p²)·n_e |
| `SectorDipole` | `Dipole` | `Length` (m)、`Angle` (°) | 偏转，并在偏转平面内产生弱聚焦和色散 |
| `CylindricalRFCavity` | `Cavity` | `Gradient` (MV/m)、`Frequency` (MHz)、`Phase` (°) | 改变粒子能量，同时有横向和纵向聚焦 |
| `Aperture` | `Aperture` | Type 为 `Circular` 时给 `Radius`；为 `Elliptical`/`Rectangular` 时给 `RadiusX`、`RadiusY` (m) | 切掉孔径外的粒子 |
| `RPLCswitch` | `RPLCswitch` | `theta`、`phi` (rad) | 切换局部坐标轴的取向 |
| `Octupole`、`QuadDoublet`、`QuadTriplet` | 参数文件不支持 | — | 只能在 Python 代码里直接创建 |

> 螺线管和 Gabor 透镜在参数文件里使用同一种强度约定（rad/m），所以 `*-Solenoid.csv` 直接沿用了 `*-Gabor.csv` 里的强度数值，方便对比两种透镜。

### 2.4 束流参数（TN §3）

在每个元件出口，对所有存活粒子计算：

- **协方差矩阵**：C₆ = ⟨φ φᵀ⟩。
- **发射度**：ε_x = √det C_x，ε_y 同理（C_x、C_y 是 2×2 子矩阵）；四维发射度 ε₄ = (det C₄)^{1/4}，六维发射度 ε₆ = (det C₆)^{1/6}。
- **Twiss 参数**：σ_u² = β_u ε_u，⟨u′²⟩ = γ_u ε_u，⟨u u′⟩ = −α_u ε_u，并满足 β_u γ_u − α_u² = 1（u = x 或 y）。

`plotBeam.py` 输出的 CSV 每行对应一个元件出口：

| 列 | 含义 |
|---|---|
| `Name` | 元件名，格式为 `设施:Stage:Section:元件:序号`，例如 `LhARA:1:Final focus:Fquad:1` |
| `nPrtcls` | 该处存活的粒子数（除以生成数就是传输效率） |
| `sigmaxy` | [σx, σy]，RMS 束斑尺寸 (m) |
| `emittance` | [εx, εy, εL（纵向）, ε4, ε6]（存活粒子少于 10 个时不计算） |
| `Twiss` | [[αx, βx, γx], [αy, βy, γy]]，β 的单位是 m |

用 Python 读取时，可以用 `ast.literal_eval` 解析这些列表（本次修复之前，`Twiss` 列在 NumPy 2 下无法解析）。

### 2.5 两种计算包络的方法

| 方法 | 脚本 / 类 | 原理 | 优点 | 局限 |
|---|---|---|---|---|
| 粒子统计 | `plotBeam.py` / `Beam` | 用追踪得到的粒子，在每个元件出口直接计算协方差矩阵 | 包含粒子丢失和色差的真实效果 | 需要足够多的粒子，统计误差约为 1/√(2N) |
| 矩阵外推 | `plotextrapolatedBeam.py` / `extrapolateBeam` | 只在起点用粒子算一次协方差矩阵，然后用 C_出 = M C_入 Mᵀ 往下游推 | 非常快，优化程序用的就是它 | 不计粒子丢失，矩阵按参考动量计算 |

**比较这两种方法的结果，可以判断线性包络近似在你的束线里是否成立。** 这本身就是一个值得写进论文的检验。

### 2.6 激光驱动源模型（TN §5）

源（`Source` 类）有 5 种模式，由参数文件里的 `SourceMode` 指定：

| SourceMode | 名称 | 参数 |
|---|---|---|
| 0 | 参数化 TNSA（LhARA 默认） | 见下表 |
| 1 | 高斯 | `SigmaX`、`SigmaY` (m)、`MinCTheta`、`MeanEnergy`、`SigmaEnergy` (MeV) |
| 2 | 平坦 | `SigmaX`、`SigmaY`、`MinCTheta`、`Kmin`、`Emax` |
| 3 | 从文件读入 | 无参数；用 `runBeamSim.py -i <文件>` 指定，BDSIM 格式再加 `-z True` |
| 4 | 均匀圆盘 | `MeanEnergy`、`SigmaEnergy`、`Radius` |

**TNSA 模型**（mode 0）的主要内容：

- **能谱**（Fuchs et al., Nature Physics 2, 48 (2006)）：dN/dK ∝ exp(−√(2K/T_e)) / √(2K T_e)，在截止能量 K_max 以上为零。
- **热电子温度 T_e**：默认由激光参数按有质动力定标（ponderomotive scaling）计算，并按 Dover et al., HEDP 37 (2020) 100847 对小焦斑做修正。
- **截止能量 K_max**：解方程
  t_laser/t₀ = X(1 + 1/(2(1−X²))) + ¼ ln((1+X)/(1−X))，其中 X = √(K_max/K_∞)。
  代码在 `Source.DurationBYt0equation` 里用 `scipy.optimize.fsolve` 求解，本次修复的 NumPy 兼容问题就出在这里。
- **角分布**：令 r = sin θ_S，r′ 按 g(r′) = 3/(4r′ₘ²)·(r′ₘ² − r′²) 抽样，分布宽度随能量变化（由 `SigmaThetaS0`、`SlopeThetaS` 控制）；方位角 φ 在 0–2π 内均匀分布。
- **横向位置**：x、y 独立的高斯分布。

TNSA 模式的参数（参数文件里不写就用默认值，默认值见 `Source.setDEFAULTparams`）：

| 参数名（CSV） | 含义 | 默认值 |
|---|---|---|
| `Wavelength` | 激光波长 | 0.8 µm |
| `Power` | 激光功率 | 2.5×10¹⁴ W（脉冲能量 10 J ÷ 脉宽 40 fs） |
| `Strehl ratio` | Strehl 比（焦斑质量） | 0.6 |
| `r0` | 焦斑半径 | 1.5 µm |
| `Duration` | 脉宽 | 40 fs |
| `Thickness` | 靶厚 | 5 µm |
| `DivAngle` | 热电子发散半角 | 25° |
| `Te` | 热电子温度 | 由激光参数计算 |
| `Kmin` / `Kmax` | 最小 / 最大动能 | 1 MeV / 由上面的方程求解 |
| `SigmaThetaS0`、`SlopeThetaS` | 角分布宽度的截距和斜率 | 20°、15° |
| `rpmax` | 最大 r′ | −9999（不限制）；LhARA 参数文件里设为 0.005 |

> 研究传输效率时要注意：`rpmax` 会限制源的角分布范围。效率的“分母”是生成的粒子数，所以先弄清 `rpmax` 怎么影响这个分母，再解读效率数字（见 TN §5 角分布部分和 `Source` 类的代码）。

---

## 3. 代码结构与数据流

### 3.1 核心模块（`01-Code/`）

| 模块 | 主要的类 / 函数 | 职责 |
|---|---|---|
| `BeamLine.py` | `BeamLine`（单例） | 读入参数文件（pandas）；解析设施（`parseFacility`）和源（`parseSource`）；逐个创建元件；推进参考粒子；`trackBeam` 生成并追踪粒子 |
| `BeamLineElement.py` | `BeamLineElement` 基类，以及 `Facility`、`Source`、`Drift`、`Aperture`、`FocusQuadrupole`、`DefocusQuadrupole`、`SectorDipole`、`Solenoid`、`GaborLens`、`CylindricalRFCavity`、`RPLCswitch`、`Octupole`、`QuadDoublet`、`QuadTriplet` | 元件的位置、方向、准直误差（`rStrt`、`vStrt`、`drStrt`、`dvStrt`）；传输矩阵；`Transport()` 追踪一个粒子并判断丢失 |
| `Particle.py` | `Particle`、`ReferenceParticle`、`UnstableParticle`，以及粒子种类 `proton`、`twelveC6`、`pion`、`muon`、`electron`、`neutrino` | 保存每个元件出口的迹空间 / 相空间；作图（`plotTraceSpaceProgression` 等）；不稳定粒子的衰变 |
| `Beam.py` | `Beam`、`extrapolateBeam` | 从数据文件计算协方差矩阵、σ、发射度、Twiss；输出 CSV 报告和图 |
| `BeamIO.py` | `BeamIO` | 二进制数据文件的读写（格式带版本号，读取时兼容旧版本，`12-Data4Tests` 里有 V1–V8 的样例），以及 BDSIM 文本格式 |
| `Simulation.py` | `Simulation`（单例） | 把“建束线 → 生成 N 个事件 → 写文件”包装起来：`Simulation(...).RunSim()` |
| `UserFramework.py` | `startAnalysis`、`handleFILES`、`EventLoop` | 用户分析框架：解析命令行，打开文件，跑事件循环并调用用户钩子 |
| `visualise.py` | `visualise` | 画束线布局和粒子轨迹（实验室系或 RPLC） |
| `PhysicalConstants.py` | `PhysicalConstants` | 物理常数和粒子性质（质量、电荷等） |
| `pionDECAY.py`、`muonDECAY.py` | — | π、μ 衰变运动学（较新的功能） |
| `Report.py`、`LaTeX.py` | — | CSV / LaTeX 表格输出 |
| `LhARALinearOptics.py` | `version()` | 读取 git 版本（tag 和 commit），写进数据文件 |

### 3.2 数据流

```
11-Parameters/*.csv
      │  pandas 读入
      ▼
BeamLine（单例）── parseFacility / parseSource / 逐行创建元件
      │            └─ ReferenceParticle 先走一遍，记录每个元件处的 p₀、位置和方向
      ▼
Simulation.RunSim() ──► BeamLine.trackBeam(N)
      │                    对每个事件：Source 抽样得到 φ₀
      │                    对每个元件：φ ← element.Transport(φ)；丢失就停止
      │                    Particle 记录每个元件出口的 φ
      ▼
BeamIO 写二进制 .dat（文件头：格式版本、代码的 git 版本、束线定义；之后逐个事件写入）
      │
      ├─► readBeamSim.py           逐元件的迹空间散点图
      ├─► plotBeam.py              Beam：σ、ε、Twiss 随 s 的变化，加 CSV
      ├─► plotextrapolatedBeam.py  extrapolateBeam：C′ = M C Mᵀ 外推
      ├─► visualiseBeam.py         束线布局和粒子轨迹
      ├─► 2BDSIM.py                导出 BDSIM 格式
      └─► UserFramework + 你自己的 UserAnal 钩子（31-UserDirectory）
```

### 3.3 写自己的代码时需要知道的几点

- **单例**：`BeamLine` 和 `Simulation` 在一个 Python 进程里只能有一个实例。要在循环里反复重建束线（比如做优化），先清理：
  `BL.BeamLine.cleaninstance()`、`BLE.BeamLineElement.cleaninstances()`、`Prtcl.Particle.cleanAllParticles()`。参考 `31-UserDirectory/01-Code/EnvelopeOptimisation.py` 里的 `setupIteration`。
- **向束线追加元件**的标准写法：创建元件 → `BL.BeamLine.addBeamLineElement(元件)` → `参考粒子.setReferenceParticle(元件)`。参考 `EnvelopeOptimisation.setBeamLine`。
- **元件位置**：下一个元件的起点 = 上一个元件的 `getrStrt() + getStrt2End()`，方向取 `getvEnd()`。
- **随机数**：源抽样用的是 Python 的 `random` 模块，`Simulation` 创建时用当前时间（秒）作种子，所以每次运行结果都有统计涨落。需要可重复的结果时，在创建 `Simulation` 之后、`RunSim()` 之前调用 `random.seed(固定整数)`。
- **可追溯性**：数据文件头里记录了写文件时代码的 git 版本。所以做正式结果前**先提交代码**，以后就能查到每个结果是哪个版本的代码算出来的。

---

## 4. 参数文件详解（11-Parameters）

### 4.1 格式（TN §6）

CSV 文件共 8 列，**从上到下按顺序**读入，每个元件占连续的若干行，每行一个参数：

| 列 | 含义 |
|---|---|
| `Stage` | 整数，阶段编号（只用于生成元件名） |
| `Section` | 分段名，例如 Capture、Arc（用于生成元件名） |
| `Element` | 元件类型：`Global`、`Source`、`Drift`、`Aperture`、`Fquad`、`Dquad`、`Gabor lens`、`Solenoid`、`Dipole`、`Cavity`、`RPLCswitch` |
| `Type` | 元件的子类型，例如 Aperture 的 `Circular`；准直误差行用 `Shift` / `Tilt` |
| `Parameter` | 参数名，例如 `Length`、`kq` |
| `Value` | 数值 |
| `Unit` | 单位，**仅作记录，代码不读这一列**。所以数值必须用代码约定的单位（见 2.3 和 2.6 节的表） |
| `Comment` | 注释 |

文件开头的设施（Facility）部分：

```
Stage,Section,Element,Type,Parameter,Value,Unit,Comment
0,Facility,Global,Name,Name,LhARA,,
0,Facility,Global,Reference particle,Kinetic energy,15,MeV,
0,Facility,Global,Reference particle,Species,12C6,,        ← 可选：粒子种类，默认是质子
0,Facility,Global,Vacuum chamber,Mother volume radius,0.035,m,   ← 真空室半径，超出即丢失
```

参考粒子的动能决定了磁刚度 Bρ，进而决定各元件的设计强度。换成离子（例如 `12C6`，即 ¹²C⁶⁺）时，记得同时改参考能量和各元件强度。

### 4.2 准直误差（做公差研究用）

在某个元件的参数行**后面紧接着**加若干行：Stage、Section、Element 都和这个元件相同，Type 写 `Shift` 或 `Tilt`。

```
1,Final focus,Fquad,,Length,0.114,m,
1,Final focus,Fquad,,kq,20.708,m-2,
1,Final focus,Fquad,Shift,dx,0.001,m,水平位移 1 mm
1,Final focus,Fquad,Shift,dy,0.0005,m,竖直位移 0.5 mm
1,Final focus,Fquad,Tilt,alphaE,0.002,rad,欧拉角（还有 betaE、gammaE）
```

这样写会设置该元件的 `drStrt`（位移）和 `dvStrt`（倾斜），已实测有效。

### 4.3 现有的参数文件

| 文件 | 内容 |
|---|---|
| `LhARABeamLine-Params-LsrDrvn-Gabor.csv` | **LhARA Stage 1 主配置**：TNSA 源 + Gabor 透镜 |
| `LhARABeamLine-Params-LsrDrvn-Solenoid.csv` | 同上，用螺线管代替 Gabor 透镜 |
| `LhARABeamLine-Params-Gauss-Gabor.csv`、`...-Gauss-Solenoid.csv` | 高斯源 |
| `LhARABeamLine-Params-Gauss-Solenoid-12C6.csv` | 碳离子 |
| `LhARA-Stage1-BeamDelivery.csv` | 从文件读入粒子，只模拟末端聚焦（当前优化工作用的配置） |
| `LIONBeamLine-Params-*.csv`、`DRACOBeamLine-Params-LsrDrvn.csv` | 已有的激光离子加速装置 LION（慕尼黑）和 DRACO（德累斯顿）的束线，用于和实验数据比较 |
| `HorizLeft/HorizRight/VerticalUp/VerticalDown-*.csv` | 不同方向的偏转，用于测试坐标变换 |
| `proton-*.csv`、`carbon-*.csv`、`decayCHAIN*.csv`、`pion4Tests.csv`、`Dummy4Tests.csv` 等 | 测试用 |
| `Template.csv` | 空模板 |

### 4.4 编辑建议

- 同名的 `.xlsx` 是原作者用 Excel 编辑的源文件，改完后另存为 CSV。也可以直接编辑 CSV。
- **不要直接改原始参数文件**：复制一份，起一个能说明改动的名字（例如 `LhARA-BD-kq22.csv`），并提交到 git。
- 想快速扫描某个参数，用 5.5 节的 `scanParameter.py`，它会自动生成临时参数文件。

---

## 5. 运行：从生成粒子到出图

### 5.1 环境配置（本地，只需做一次）

```bash
git clone https://github.com/cid02394857/LhARAlinearOptics.git
cd LhARAlinearOptics
git checkout claude/busy-euler-mjur5i
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

每次打开新终端：

```bash
cd LhARAlinearOptics
source venv/bin/activate
source startup.bash        # 必须在 bash 里、在仓库根目录执行
```

macOS 默认的 shell 是 zsh，先输入 `bash` 切换。Windows 推荐用 WSL；如果用 Git Bash，先执行 `unset HOMEPATH`（见第 7 节）。

### 5.2 流程 A：从激光源开始模拟整条 Stage 1 束线

```bash
# 1) 生成并追踪 1000 个粒子（约 5 秒）
python 03-Scripts/runBeamSim.py -b 11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv \
       -o 99-Scratch/LhARA.dat -n 1000

# 2) 计算束流包络并写 CSV（约 5 秒）→ 99-Scratch/BeamProgressionPlot.pdf
python 03-Scripts/plotBeam.py -i 99-Scratch/LhARA.dat -n 1000 -o 99-Scratch/LhARA-beam.csv

# 3) 用矩阵外推包络（约 5 秒）；同样写 BeamProgressionPlot.pdf（会覆盖第 2 步的图）
python 03-Scripts/plotextrapolatedBeam.py -i 99-Scratch/LhARA.dat -n 1000 -o 99-Scratch/LhARA-exbeam.csv

# 4) 逐元件的迹空间散点图（约 3 分钟）
#    → 99-Scratch/ParticleProgressionPlot.pdf、ParticleLongiProgressionPlot.pdf
python 03-Scripts/readBeamSim.py -i 99-Scratch/LhARA.dat -n 1000

# 5) 束线布局和粒子轨迹
python 03-Scripts/visualiseBeam.py -i 99-Scratch/LhARA.dat -o 99-Scratch/LhARA-vis.pdf -n 100

# 6) 导出 BDSIM 格式（注意：-l 选项目前被解析了但没有被使用，见 PROGRESS.md 待办）
python 03-Scripts/2BDSIM.py -i 99-Scratch/LhARA.dat -o 99-Scratch/LhARA-bdsim.dat -n 100
```

要在不同参数之间比较 `plotBeam.py` 的图，每次运行后把 `BeamProgressionPlot.pdf` 改个名字保存下来，否则会被覆盖。

### 5.3 流程 B：束流输运段（从文件读入粒子）

```bash
python 03-Scripts/runBeamSim.py -b 11-Parameters/LhARA-Stage1-BeamDelivery.csv \
       -i 12-Data4Tests/BDSIM-Data4Tests.dat -z True -o 99-Scratch/BD.dat -n 1000
python 03-Scripts/plotBeam.py -i 99-Scratch/BD.dat -n 1000 -o 99-Scratch/BD-beam.csv
```

`-z True` 表示输入文件是 BDSIM 格式（每行是一个粒子的坐标，用空格分隔）。输入文件也可以是本代码写出的 `.dat` 文件，这时不加 `-z`。

### 5.4 流程 C：用户分析框架（`31-UserDirectory`）

用户目录的设计目的是把你自己的分析代码和主仓库的代码分开：

```bash
cp -r 31-UserDirectory ~/LhARA-work        # 拷贝到自己的工作目录
cd ~/LhARA-work
source <仓库路径>/venv/bin/activate
source user_startup.bash -p <仓库路径>      # 会自动调用主仓库的 startup.bash
$LhARAOpticsPATH/03-Scripts/runBeamSim.py -b 11-Parameters/LIONBeamLine-Params-Gauss.csv -o 99-Scratch/test.dat -n 1000
03-Scripts/visualiseBeam.py -i 99-Scratch/test.dat -o 99-Scratch/test.pdf -n 100
03-Scripts/UserAnalysis.py -i 99-Scratch/test.dat -n 100
```

- `01-Code/UserAnal.py`：分析模板，提供三个钩子：`UserInit`（初始化）、`UserEvent(粒子)`（事件循环里对每个粒子调用）、`UserEnd`（最后输出统计和图）。你的分析就写在这里。
- `01-Code/EnvelopeOptimisation.py` + `03-Scripts/EnvelopeAnalysis.py`：**四极铁俘获段的随机搜索优化**（最近一次提交的内容，约 10 秒）：

  ```bash
  03-Scripts/EnvelopeAnalysis.py -i 21-DataSets/PoPLaRsource.dat -o 99-Scratch/EO.dat -n 1000
  ```

  它的逻辑：
  1. 读入源粒子分布，在后面接一个“漂移 → 聚焦四极铁 → 漂移 → 散焦四极铁 → 漂移到 z = 2 m”的俘获段（PoPLaR 束线的配置）。
  2. **代价函数** = π σx σy × (σx/σy 和 σy/σx 中较大的那个)，即在 z = 2 m 处束斑面积尽量小，同时尽量圆。
  3. 迭代 2000 次：每次在当前最优参数附近加高斯随机扰动（漂移长度、四极铁长度和梯度），约束是前四个元件总长 < 1 m、梯度 ≥ 90 T/m；代价变小就接受。这是最简单的随机搜索，相当于没有步长自适应的 (1+1) 进化策略。
  4. 每次迭代用 `extrapolateBeam`（矩阵外推）计算包络，所以很快。最终的包络图写到 `99-Scratch/BeamProgressStrt.pdf`（优化前）和 `BeamProgressEnd.pdf`（优化后）。

### 5.5 参数扫描示例：`scanParameter.py`

这是本次新增的示例脚本，在用户目录 `03-Scripts/` 下。它每次修改参数文件里的一个参数，依次运行 `runBeamSim.py` 和 `plotBeam.py`，把束线末端的粒子数和束斑尺寸汇总起来：

```bash
# 在用户工作目录（或仓库根目录）下：
03-Scripts/scanParameter.py -b 11-Parameters/LhARA-Stage1-BeamDelivery.csv \
    -e Fquad -p kq -v 18,19,20,21,22,23,24 -n 1000 \
    -i 12-Data4Tests/BDSIM-Data4Tests.dat -z
# 输出：99-Scratch/scan-Fquad-kq.csv 和 scan-Fquad-kq.pdf
```

- `-e`、`-p` 对应参数文件的 Element 列和 Parameter 列；有多行匹配时，用 `-k` 指定第几个（默认第 1 个）。
- 文件依次在当前目录、`$HOMEPATH`、`$LhARAOpticsPATH` 下查找，所以主仓库的参数文件和数据文件可以直接用相对路径。
- 实测（每点 500 个事件）：kq = 18、20、22、24 m⁻² 时，末端 σx = 7.8、4.5、1.4、2.3 mm，可见水平束斑在 kq ≈ 22 m⁻² 附近最小（参数文件里的标称值是 20.708）。

### 5.6 跑全部测试

```bash
bash 02-Tests/runAllTests.bash                 # 约 20 分钟；日志在 99-Scratch/TestLogs/
bash 02-Tests/runAllTests.bash DriftTst        # 只跑指定的测试
```

**改了 `01-Code` 里的任何代码之后，都应该跑一遍。** 元件测试会把传输结果和预先算好的值比较，这些值可以用 `21-Spreads4Tests/` 里的电子表格核对。

---

## 6. 如何用这套代码做有效的科学研究

### 6.1 第一步：复现与验证（建议第 1–2 周）

1. 跑通 `runAllTests.bash`，确认环境正确。
2. 用流程 A 跑 LhARA 主配置，把结果和 `41-ReferencePlots/`（2024-02 的参考图）比较。数量级、包络形状应该一致；统计涨落会带来细节上的差别。
3. 读 TN 第 4 章，挑一两个元件（比如四极铁），对照 `21-Spreads4Tests/` 里的电子表格自己手算一次传输矩阵，确认你理解了约定（单位、符号、Bρ 的计算）。
4. 用 2.5 节的两种包络方法分别计算同一条束线，比较两者的差别，弄清线性近似在哪里失效（通常在源附近和能量选择段）。

### 6.2 可以做的研究方向

| 方向 | 研究的问题 | 怎么做 | 主要用到的文件 |
|---|---|---|---|
| **A. 终端站束流输运优化**（当前重点） | 怎样的末端聚焦参数能在终端站得到足够大、足够均匀的照射野，同时保证传输效率？ | 先用 `scanParameter.py` 扫单个参数，弄清楚趋势；再改进 `EnvelopeOptimisation.py`：把代价函数换成“束斑尺寸接近目标值 + 均匀度 + 传输效率”，把随机搜索换成 `scipy.optimize`（Nelder–Mead、差分进化等），比较收敛速度和结果 | `LhARA-Stage1-BeamDelivery.csv`、`EnvelopeOptimisation.py`、`scanParameter.py` |
| **B. 源参数敏感性** | 激光功率、焦斑、脉宽、靶厚变化时，终端站的束流（能谱、束斑、粒子数）怎么变？ | 在 TNSA 源参数里加 `Power`、`r0`、`Duration`、`Thickness` 等行，逐个扫描；同时记录 K_max 和能谱 | 主配置 CSV、`ExpSourceTst.py`（源分布的诊断图） |
| **C. 能量选择性能** | 准直器半径、RF 腔相位和梯度，怎样影响终端站的能散和传输效率？ | 扫 `Aperture` 的 `Radius`、`Cavity` 的 `Phase`/`Gradient`；用 `readBeamSim.py` 的纵向图看能散 | 主配置 CSV |
| **D. 公差研究** | 元件准直误差、强度误差有多大时，终端站的束流会超出指标？ | 用 4.2 节的 Shift/Tilt 行，对每个元件加随机误差，做多次蒙特卡洛（每次换随机种子），统计束斑中心的偏移和尺寸的分布 | 主配置 CSV |
| **E. Gabor 透镜与螺线管对比** | 两种透镜在相同强度约定下，发射度增长、x–y 耦合、色差有什么区别？ | 用 `*-Gabor.csv` 和 `*-Solenoid.csv` 跑同样的事件数，比较 εx、εy、ε4 和束斑 | 两组参数文件 |
| **F. 离子束** | 碳离子（¹²C⁶⁺）和质子相比，需要怎样调整束线？ | 用 `*-12C6.csv`，注意磁刚度的变化 | `LhARABeamLine-Params-Gauss-Solenoid-12C6.csv`、`*-12C6-Tst.py` |
| **G. 与 BDSIM 交叉验证** | 线性光学的结果和全追踪差多少？ | 用 `2BDSIM.py` 把同一批粒子导给 BDSIM，比较终端站的分布 | `2BDSIM.py`、`-z` 选项 |

### 6.3 科研规范（让结果站得住脚）

- **统计误差**：σ 的相对统计误差约为 1/√(2N)。比较两个方案时，差别要明显大于统计误差。正式结果用 1 万个以上的事件，并检查事件数加倍时结论是否不变（收敛性）。
- **可重复性**：固定随机种子（3.3 节）；参数文件和代码都提交到 git，数据文件会自动记录代码版本。
- **一次只改一个量**：扫描时只改一个参数，其他保持在基准值，并且始终保留一个“基准配置”的结果作对照。
- **交代近似**：写论文时说明 1.4 节列出的近似，并用 2.5 节的对比和 BDSIM 交叉验证说明这些近似在你的问题里是否够用。
- **做记录**：像原作者那样在 `91-Notes/` 写日期日志（做了什么、命令、结果、结论），对写论文帮助很大。
- **改代码之后跑测试**：`bash 02-Tests/runAllTests.bash`。

---

## 7. 常见问题

| 现象 | 原因和解决办法 |
|---|---|
| `source startup.bash` 报 “environment appears to already be configured” | 当前 shell 里已经 source 过一次（`HOMEPATH` 已设置）。新开一个终端，或者先 `unset HOMEPATH`。**Windows 的 Git Bash 自带系统变量 `HOMEPATH`**，必须先 unset；推荐用 WSL |
| `ModuleNotFoundError` | 没有激活 venv，或者没有 source `startup.bash`（它负责把 `01-Code` 加进 `PYTHONPATH`） |
| 找不到 `11-Parameters/...` 或 `99-Scratch` | 测试脚本要在仓库根目录下运行；用户脚本要在你的工作目录下运行 |
| 图形窗口弹出、程序卡住 | 执行 `export MPLBACKEND=Agg`，只保存图片、不弹窗口（`runAllTests.bash` 已自动设置） |
| 长选项（如 `--bfile`）不起作用 | 部分脚本的长选项定义不完整，请用短选项（`-b`、`-i`、`-o`、`-n`） |
| 每次运行结果不一样 | 随机种子取自当前时间，见 3.3 节 |
| `readBeamSim.py` 很慢 | 它要为每个元件画散点图；作图时减小 `-n` |
| 没装 LaTeX | 不影响使用；装了 LaTeX 的话，图里的公式会用 LaTeX 渲染 |
| 依赖版本 | 2026-10 用最新版验证过：numpy 2.5、pandas 3.0、scipy 1.18、matplotlib 3.11 |

---

## 8. 附录：文件速查与 TN 章节对照

### 8.1 你最常用的文件

| 目的 | 文件 |
|---|---|
| 改束线 | `11-Parameters/*.csv` |
| 跑模拟 | `03-Scripts/runBeamSim.py` |
| 看包络 / 传输效率 | `03-Scripts/plotBeam.py`、`plotextrapolatedBeam.py` |
| 看粒子分布 | `03-Scripts/readBeamSim.py`、`visualiseBeam.py` |
| 写自己的分析 | `31-UserDirectory/01-Code/UserAnal.py` |
| 优化 | `31-UserDirectory/01-Code/EnvelopeOptimisation.py`、`03-Scripts/EnvelopeAnalysis.py` |
| 参数扫描 | `31-UserDirectory/03-Scripts/scanParameter.py` |
| 验证 | `02-Tests/runAllTests.bash`、`21-Spreads4Tests/`、`41-ReferencePlots/` |

### 8.2 TN 章节对照

| TN 章节 | 内容 | 对应代码 |
|---|---|---|
| §2 坐标系 | 实验室系、RPLC | `BeamLineElement`（`Tilt2RPLC` 等）、`RPLCswitch` |
| §3 相空间与迹空间 | φ 的定义、协方差、发射度、Twiss | `Particle`、`Beam` |
| §4 传输矩阵 | 各元件的矩阵 | `BeamLineElement` 的各子类（`setTransferMatrix`） |
| §5 源 | TNSA 模型 | `BeamLineElement.Source`、`BeamLine.parseSource` |
| §6 束线描述与脚本 | CSV 格式、`03-Scripts` | `BeamLine`、`03-Scripts/*.py` |
| 附录 A1 | 各模块、类和数据格式的说明 | `01-Code/*.py` 文件开头的说明文档 |
| 附录 A2 | 安装与运行 | `README.md`、`startup.bash` |
