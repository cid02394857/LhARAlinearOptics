# CLAUDE.md — LhARAlinearOptics 项目记忆

> 本文件供 Claude Code（以及其他 AI 编程助手）在每次新对话开始时自动读取。
> **与用户交流请使用中文。** 用户是硕士生，项目是 LhARA 第一阶段（Stage 1）束线的代码模拟。
> 开始工作前请先读 `PROGRESS.md`：里面有进度日志、已修复问题和待办事项。

## 项目概况

- LhARA = Laser-hybrid Accelerator for Radiobiological Applications（用于放射生物学研究的激光混合加速器）。
- 本仓库是线性束流光学（linear beam optics）模拟代码：用传输矩阵（transfer matrix）在迹空间（trace space）中，
  把激光驱动源产生的质子/离子逐个元件地追踪到束线末端，并统计束流包络、发射度、Twiss 参数等。
- 上游仓库：`ImperialCollegeLondon/LhARAlinearOptics`；用户 fork：`cid02394857/LhARAlinearOptics`（**没有 git tag**）。
- 方法论文档：`00-Documentation/01-WriteUp/`（LaTeX 源码，编译版就是 `LhARA-TN-2024-07.pdf`）。
- 中文代码解读与科研使用指南：`00-Documentation/CodeGuide-CN.md`。

## 目录结构

| 目录 | 内容 |
|---|---|
| `01-Code/` | 核心类库：`BeamLine`、`BeamLineElement`（全部元件类）、`Particle`、`Beam`、`BeamIO`、`Simulation`、`UserFramework` 等 |
| `02-Tests/` | 测试脚本，基本上每个类一个（作者的标准测试清单见 `91-Notes/2026-09-02.txt`） |
| `03-Scripts/` | 面向用户的命令行脚本：`runBeamSim.py`、`readBeamSim.py`、`plotBeam.py`、`plotextrapolatedBeam.py`、`visualiseBeam.py`、`2BDSIM.py` |
| `11-Parameters/` | 束线参数 CSV（同名 xlsx 是编辑用源文件） |
| `12-Data4Tests/` | 测试用数据文件（二进制 BeamIO 格式，外加 BDSIM 文本格式） |
| `31-UserDirectory/` | 用户分析框架示例，包括四极铁优化 `01-Code/EnvelopeOptimisation.py` |
| `91-Notes/` | 原作者的开发日志 |
| `99-Scratch/` | 输出目录（已 gitignore，由 `startup.bash` 自动创建） |

## 环境配置（本地，bash）

```bash
python3 -m venv venv && source venv/bin/activate
pip install -r requirements.txt
source startup.bash          # 必须在仓库根目录、在 bash 里执行
export MPLBACKEND=Agg        # 可选：没有图形界面时用
python 02-Tests/BeamLineTst.py
```

- `startup.bash` 会设置 `HOMEPATH`、`LhARAOpticsPATH`、`PYTHONPATH`（加入 `01-Code`）和 `REPORTPATH`（`99-Scratch`）。
  如果 `HOMEPATH` 已经存在，脚本会报错并退出：同一个 shell 里重复 source 会这样；**Windows Git Bash 的 `HOMEPATH` 是系统变量**，需要先 `unset HOMEPATH`（推荐用 WSL）。
- 测试脚本要在**仓库根目录**下运行（脚本里用的是相对路径 `11-Parameters/...`、`12-Data4Tests/...`）。
- 依赖：GitPython、matplotlib、numpy、pandas、pylorentz、scipy（2026-10 已验证可用的最新版本：numpy 2.5.3、pandas 3.0.5、scipy 1.18.1、matplotlib 3.11.2）。
- LaTeX 不是必需的：没装 LaTeX 时，作图会自动退回 matplotlib 自带的公式渲染。

## 修改代码时的约定

- 跟随现有代码风格：类方法 `getX` / `setX`、调试输出放在 `if self.getDebug():` 块里、注释用 `#..` 开头、续行用 `\`。
- 二进制数据格式（`BeamIO`）是带版本号的（目前是 V8），读取时要兼容旧版本；改格式时要同时改读写两端，并更新 `12-Data4Tests` 和 `BeamIOTst-rd.py`。
- 改完代码后，跑一遍完整测试（见 `PROGRESS.md` 里的“如何跑全部测试”）。
- 提交信息用英文，简洁说明改了什么、为什么改。
