# HANDOFF.md — 仓库文件与代码可读性检查：交接文档

> 编写日期：2026-10-28。检查和修复在 2026-10-23 至 2026-10-26 之间完成。
> 接手前请先读 `CLAUDE.md`（项目记忆、环境配置）和 `PROGRESS.md`（进度日志）。本文件只总结这次检查。
> 行号均为当前分支 `claude/busy-euler-mjur5i` 上的行号；已修改过的位置，括号里注明起点提交中的原行号。

---

## 1. 检查范围

**仓库与版本**

| 项 | 内容 |
|---|---|
| 仓库 | `cid02394857/LhARAlinearOptics`（`ImperialCollegeLondon/LhARAlinearOptics` 的 fork，**没有 git tag**） |
| 分支 | `claude/busy-euler-mjur5i`（尚未合并进 `main`） |
| 检查起点 | `d75da9c` “Fixes to allow analysis of quad optimisation for LhARA beam delivery, debug printout.”（检查开始时等于 `origin/main`） |
| 检查结束时 | `7831e16`（本文件在其后单独提交） |
| 运行环境 | Linux 云容器，Python 3.13.16；numpy 2.5.3、pandas 3.0.5、scipy 1.18.1、matplotlib 3.11.2、GitPython 3.2.1、pylorentz；没有安装 LaTeX。另外用 numpy 2.4.6 和 2.3.5 交叉验证过 NumPy 兼容问题 |

**看过的目录和文件**

| 范围 | 做了什么 |
|---|---|
| 全部 359 个文件 | 检查读权限、文件类型（`file`）、是否为 Git LFS 占位文件 |
| 全部 101 个 `.py`（含 `ZZ-Archive`） | 编译检查，并把 SyntaxWarning 当作错误再检查一遍 |
| `01-Code/` | 阅读 `BeamLine.py`（解析参数、`parseSource`）、`BeamLineElement.py`（`Transport`、`Source`、`Drift`、`GaborLens`、`Octupole`、`RPLCswitch` 等）、`Particle.py`、`Beam.py`、`BeamIO.py`、`Simulation.py`、`UserFramework.py`、`LhARALinearOptics.py`、`LaTeX.py` 的相关部分 |
| `02-Tests/` | 56 个测试脚本全部运行过（`ZZ-Archive` 除外） |
| `03-Scripts/` | 6 个脚本全部阅读，并用实际参数运行 |
| `31-UserDirectory/` | 按 README 拷贝到其他目录后运行 4 个脚本；阅读 README、`user_startup.bash`、`EnvelopeOptimisation.py` |
| `11-Parameters/` | 主要看了 `LhARABeamLine-Params-LsrDrvn-Gabor.csv`、`LhARA-Stage1-BeamDelivery.csv`、Solenoid 版和 12C6 版 |
| `12-Data4Tests/` | 通过 `BeamIOTst-rd.py`、`BeamIOTst-rd-BDSIM.py` 实际读取 V1–V8 二进制文件和 BDSIM 文件 |
| `00-Documentation/01-WriteUp/` | 阅读 LaTeX 源码第 1–6 章和附录 A2；确认用户上传的 `LhARA-TN-2024-07.pdf` 与仓库里的副本 md5 一致 |
| `91-Notes/` | 阅读 `2026/2025-09-05.txt`、`2026/2025-11-15.txt`、`2026-01-05.txt`、`2026-01-12.txt`、`2026-09-02.txt` |
| 根目录 | `README.md`、`requirements.txt`、`startup.bash`、`.gitignore` |
| GitHub 远程 | 用 `git fetch`/`ls-remote` 和 GitHub 接口确认可访问 |

**没有详细检查**：`01-Code/ZZ-Archive`、`02-Tests/ZZ-Archive`（只做了编译检查）；xlsx 文件内容；`21-Spreads4Tests/`；`41-ReferencePlots/`（只检查了文件类型）；附录 `A1-MCD.tex`；`91-Notes/` 中的其他日志；`92-Literature/`（只确认了 Mora 那篇 PRL）。

**整体结论**：359 个文件全部可读，没有 LFS 占位文件；101 个 `.py` 都能编译，没有语法错误或警告；PDF、xlsx、图片、数据文件的格式都正常。问题集中在**运行时**，见第 2 节。

---

## 2. 主要发现（按严重程度）

严重程度：**严重** = 程序崩溃、无法运行；**高** = 在常见环境下失败，或输出不可用；**中** = 能运行，但影响结果解读或物理正确性；**低** = 易用性、文档或潜在问题。

### 严重（3 项，均已修复）

| # | 问题 | 位置 | 状态 |
|---|---|---|---|
| S1 | `requirements.txt` 漏了 `pylorentz`，而核心模块 import 了它；只按 `requirements.txt` 安装时，几乎所有测试一启动就报 `ModuleNotFoundError` | `01-Code/Particle.py:150`（import）；`requirements.txt:5`（补上） | 已修复 |
| S2 | fork 没有 git tag，`version()` 里的 `tags[-1]` 抛出 `IndexError`；所有写数据文件的程序都经 `BeamIO` 调用它，因此全部崩溃（模拟、写束流文件等） | `01-Code/LhARALinearOptics.py:42–48`（原 39、41 行）；调用处 `01-Code/BeamIO.py:493` | 已修复 |
| S3 | NumPy ≥ 2.4 时，`math.log` 作用在 `fsolve` 传入的单元素数组上，抛出 `TypeError`；激光驱动源（LsrDrvn）配置连束线都建不起来 | `01-Code/BeamLineElement.py:6071`（原 6068）；调用处 `01-Code/BeamLine.py:664` | 已修复 |

### 高（5 项，均已修复）

| # | 问题 | 位置 | 状态 |
|---|---|---|---|
| H1 | 输出目录 `99-Scratch` 在新克隆的仓库里不存在，所有往里写文件的测试都失败 | `startup.bash:13`（新增 `mkdir -p`） | 已修复 |
| H2 | 3 个测试脚本强制 `usetex=True`，没装 LaTeX 时作图报 `RuntimeError` | `02-Tests/decayCHAIN.py:30`、`02-Tests/ExpSourceTst.py:30`、`02-Tests/ExpSourceTst_rms_KE.py:27` | 已修复（装了 LaTeX 才启用） |
| H3 | NumPy 2 下，`plotBeam.py`/`plotextrapolatedBeam.py` 输出 CSV 的 `Twiss` 列被写成 `np.float64(...)`，无法解析 | `01-Code/Beam.py:536–547` | 已修复 |
| H4 | `OctupoleTst.py` 一直失败（原作者笔记里也记为失败）：测试用的是旧版 `vStrt` 格式；`Octupole` 类没有设置 `Strt2End`/`Rot2LbEnd`，放不进束线；预期值本身前后矛盾 | `01-Code/BeamLineElement.py:3709–3710`（补的两行）；`02-Tests/OctupoleTst.py:37、39、52、104` | 已修复 |
| H5 | 脚本运行成功时退出码为 1（用法错误时反而是 0），批量运行或扫描时无法判断成败 | 各文件最后一行：`03-Scripts/2BDSIM.py:117`、`plotBeam.py:113`、`plotextrapolatedBeam.py:109`、`readBeamSim.py:132`、`runBeamSim.py:127`、`visualiseBeam.py:176`；`31-UserDirectory/03-Scripts/EnvelopeAnalysis.py:104`、`UserAnalysis.py:46`、`visualiseBeam.py:219`；`02-Tests/decayCHAIN-vis.py:214` | 已修复 |

### 中（6 项，未修改：涉及物理模型，需要和导师或原作者确认）

| # | 问题 | 位置 | 状态 |
|---|---|---|---|
| M1 | 每个元件入口都有纵向截断 \|z\| > 5 m。主配置 5 万事件中，**91.4%** 的“损失”来自它，而不是物理孔径；全束线传输率和末端能谱主要由它决定 | `01-Code/BeamLineElement.py:635`（`BeamLineElement.Transport`）、`:5252`（`CylindricalRFCavity.Transport`） | 未修改 |
| M2 | 偏转段入口的 `RPLCswitch` 丢弃了 3.5% 的生成粒子，而且都是 13–17 MeV 的“好”粒子，原因未查明 | `01-Code/BeamLineElement.py:7753`（`RPLCswitch.Transport`）；参数见 `11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv:61` 起 | 未修改 |
| M3 | 线性输运中，远离参考动量的粒子横向位置被缩放了 p/p₀，原因有三：源处 x′ = xₚ·p/p₀；漂移段用参考动量下的固定矩阵；Gabor 透镜用按粒子动量计算的物理斜率矩阵。比值类的量（x′/x、焦点位置）不受影响；5 MeV 的束斑偏小到 0.57 倍，30 MeV 偏大到 1.42 倍 | `01-Code/BeamLineElement.py:6241`（源处 x′）、`:638`（只有这里列出的元件按粒子动量重算矩阵）、`:4513`（Gabor 透镜矩阵） | 未修改 |
| M4 | 源模型先按全角度能谱抽能量，再强制 r′ < `rpmax`，不随能量调整接受率。主配置设 `rpmax = 0.005`，只相当于 TNSA 锥内约 0.08% 的质子，而且锥内能谱偏软 | `01-Code/BeamLineElement.py:5934`（抽能量）、`:5972`（`rpmax` 判断）；`11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv:6` | 未修改 |
| M5 | 能量转换效率的定标 1.2×10⁻¹⁵·I^0.75 在 2.1×10²¹ W/cm² 下得到 11.9，代码截断为 0.5；TN 注明它只适用于 I < 3.1×10¹⁹ W/cm² | `01-Code/BeamLineElement.py:5871`；`00-Documentation/01-WriteUp/05-Source/05-Source.tex:56`（行内注释） | 未修改 |
| M6 | 质子源横向尺寸取激光焦斑 r0 = 1.5 µm，而模型中的电子鞘层半径是 3.8 µm，源尺寸可能被低估 | `01-Code/BeamLineElement.py:5939` | 未修改 |

### 低（7 项）

| # | 问题 | 位置 | 状态 |
|---|---|---|---|
| L1 | Gabor 透镜强度的单位不一致：参数文件写 rad/m，调试输出写 m⁻²，代码内部实际存的是等效磁场 B₀（= Strength × Bρ，单位 T） | `11-Parameters/LhARABeamLine-Params-LsrDrvn-Gabor.csv:13`；`01-Code/BeamLineElement.py:4280` | 未修改 |
| L2 | 元件矩阵用参考粒子**最后一个记录位置**的动量计算（`iPrev = len(...) - 1`）。追踪时取到的是束线末端的参考动量，目前参考动量沿线不变，所以没有影响；一旦参考粒子被加速就会出错 | `01-Code/BeamLineElement.py:4459`（`GaborLens.setTransferMatrix`）；同样写法见 2396、2424（`FocusQuadrupole`）、2973、3001（`DefocusQuadrupole`）、3528、3557（`SectorDipole`）、4075、4105（`Solenoid`） | 未修改 |
| L3 | `2BDSIM.py` 解析了 `-l` 选项但没有使用，导出位置是写死的 | `03-Scripts/2BDSIM.py:24、39`（解析）、`:98`（写死） | 未修改 |
| L4 | getopt 长选项定义缺少 `=`，`--bfile x` 这种写法不起作用；`-z` 用 `bool(arg)` 解析，`-z False` 也被当作 True | `03-Scripts/runBeamSim.py:21`、`:47`（其他脚本有同样的长选项问题） | 未修改 |
| L5 | `startup.bash` 发现 `HOMEPATH` 已存在就退出；Windows Git Bash 里 `HOMEPATH` 是系统变量，所以直接用会失败 | `startup.bash:3` | 文档中已说明，未改代码 |
| L6 | 文档错误或过时：README 里的安装指南路径写错；用户目录 README 写成了 `startup.bash`；TN 说代码“只处理质子”；TN 写的输出文件名 `extrapolatedBeamProgress.pdf` 已过时，实际两个包络脚本都写 `BeamProgressionPlot.pdf`，后运行的会覆盖先运行的 | `README.md:7`、`31-UserDirectory/README.md:22`（已修复）；`00-Documentation/01-WriteUp/01-Intro/01-Intro.tex:9`、`06-BLspecScrpts/06-BLspecScrpts.tex:124`、`01-Code/Beam.py:908`（未修改） | 部分修复 |
| L7 | `LaTexTst.py` 按默认路径在仓库根目录写出 `LaTex.tex`，每次跑测试都会在 `git status` 里出现 | `01-Code/LaTeX.py:7`；`.gitignore:6`（新增忽略） | 已处理 |

---

## 3. 已做的修改

共 7 次提交（`ba7fa3e`、`aa46f6c`、`733d3c3`、`9831d85`、`4271392`、`655fe74`、`7831e16`），都已推送到 `origin/claude/busy-euler-mjur5i`。

**修复（对应第 2 节）**

| 文件 | 改动 | 原因 |
|---|---|---|
| `requirements.txt` | 加入 `pylorentz` | S1 |
| `01-Code/LhARALinearOptics.py` | 没有 tag 时返回 `"No tag"`；找不到 git 仓库时按文档约定的格式返回占位值 | S2 |
| `01-Code/BeamLineElement.py` | 第 6071 行 `mth.log` 改为 `np.log`；`Octupole.__init__` 补上 `setStrt2End`、`setRot2LbEnd` 两行 | S3、H4 |
| `01-Code/Beam.py` | Twiss 参数转为 Python float | H3 |
| `startup.bash` | 自动创建 `99-Scratch` | H1 |
| `02-Tests/decayCHAIN.py`、`ExpSourceTst.py`、`ExpSourceTst_rms_KE.py` | 装了 LaTeX 才启用 `usetex` | H2 |
| `02-Tests/OctupoleTst.py` | 改用新版 `vStrt` 格式，先建束线、登记参考粒子，按漂移段矩阵重算预期值 | H4 |
| `03-Scripts/` 6 个脚本、`31-UserDirectory/03-Scripts/` 3 个脚本、`02-Tests/decayCHAIN-vis.py` | 最后一行 `sys.exit(1)` 改为 `sys.exit(0)` | H5 |
| `README.md`、`31-UserDirectory/README.md` | 改正安装指南路径和脚本名 | L6 |
| `.gitignore` | 忽略 `LaTex.tex` | L7 |

**新增**

| 文件 | 用途 |
|---|---|
| `CLAUDE.md` | 项目记忆：本地 Claude Code 每次新对话会自动读取 |
| `PROGRESS.md` | 进度日志、已修复问题、待办事项 |
| `HANDOFF.md` | 本文件 |
| `02-Tests/runAllTests.bash` | 按依赖顺序一键跑全部测试，汇总 PASS/FAIL |
| `31-UserDirectory/03-Scripts/scanParameter.py` | 参数扫描示例 |
| `00-Documentation/CodeGuide-CN.md` | 中文代码解读与科研使用指南 |
| `51-Reports/2026-10-26-Source-GaborLens1/` | 主配置 5 万事件模拟报告（中英文）、分析脚本、图和结果表 |

**验证结果**
- `bash 02-Tests/runAllTests.bash`：61/61 全部通过（顺序执行约 18 分钟）。
- `03-Scripts` 的 6 个脚本、从 BDSIM 文件读入、用户目录的 4 个脚本，共 12 项检查全部通过。
- 修复 S3 前后对比：NumPy 2.3.5 只给警告，2.4.6 和 2.5.3 报错；修复后 2.5.3 下通过。

---

## 4. 未完成事项与建议的下一步

**未完成**
1. 用户还没在自己的电脑上配置环境、跑通测试。步骤见 `00-Documentation/CodeGuide-CN.md` 第 5.1、5.6 节。
2. 第 2 节“中”级的 M1–M6 都没有修改。它们涉及物理模型的设计取舍，需要先确认原作者的本意。
3. 低级问题 L1–L4 没有修改，原因同上（比如 `2BDSIM.py` 的 `-l` 本意是导出哪个位置）。
4. 分支 `claude/busy-euler-mjur5i` 还没有合并进 `main`。新开的云端会话默认在 `main` 上，看不到这些修改。

**建议的下一步**（按优先级）
1. **在本地跑通**：`git checkout claude/busy-euler-mjur5i` → 建 venv → `pip install -r requirements.txt` → `source startup.bash` → `bash 02-Tests/runAllTests.bash`。
2. **合并分支**：确认无误后合并进 `main`（在 GitHub 上开 PR 合并，或在本地 merge 后 push）。
3. **和导师确认三件事**：M1 纵向截断的物理含义（不确认的话不要用全束线传输率下结论）；M4 `rpmax` 截断是否符合研究目的；LhARA 实际使用的激光参数（目前用的是代码默认值）。
4. **查 M2**：从 `01-Code/BeamLineElement.py:7753`（`RPLCswitch.Transport`）入手，看 13–17 MeV 的粒子为什么在这里被丢弃。
5. **向上游反馈**：S1–S3、H1–H5 是通用的崩溃或兼容性修复，可以整理后向 `ImperialCollegeLondon/LhARAlinearOptics` 提 PR。
6. **科研**：见 `51-Reports/2026-10-26-Source-GaborLens1/report-CN.md` 第 7 节和 `00-Documentation/CodeGuide-CN.md` 第 6 节。
