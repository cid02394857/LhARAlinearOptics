# PROGRESS.md — 进度日志与待办事项

> 新对话（本地或云端）请先读 `CLAUDE.md`，再读本文件。每完成一项工作就更新本文件并提交。

## 如何在本地接上进度

```bash
git clone https://github.com/cid02394857/LhARAlinearOptics.git   # 已有本地仓库的话跳过这步
cd LhARAlinearOptics
git fetch origin
git checkout claude/busy-euler-mjur5i     # 本次云端会话的全部工作都在这个分支上
```

然后按 `CLAUDE.md` 的“环境配置”装好环境，在仓库根目录启动 Claude Code，对它说：
“读一下 PROGRESS.md，接着做待办事项。”

## 已修复的问题（2026-10-23）

| # | 问题 | 位置 | 修复方法 |
|---|---|---|---|
| 1 | `requirements.txt` 漏了 `pylorentz`，核心模块一 import 就报错 | `requirements.txt` | 补上 `pylorentz` |
| 2 | NumPy ≥ 2.4 下，激光驱动源（LsrDrvn）配置报 `TypeError`：对 `fsolve` 传进来的数组用了 `math.log` | `01-Code/BeamLineElement.py:6068` | `mth.log` 改成 `np.log` |
| 3 | 仓库没有 git tag（fork 就是这样）时，所有写数据文件的程序都会 `IndexError` 崩溃 | `01-Code/LhARALinearOptics.py` `version()` | 没有 tag 时返回 `"No tag"`；找不到 git 仓库时按文档约定的格式返回占位值 |
| 4 | 输出目录 `99-Scratch` 在新克隆的仓库里不存在，写输出的测试会失败 | `startup.bash` | 自动 `mkdir -p "${REPORTPATH}"` |
| 5 | 没装 LaTeX 时，`decayCHAIN.py`、`ExpSourceTst.py`、`ExpSourceTst_rms_KE.py` 作图报 `RuntimeError` | 这三个测试脚本 | 装了 LaTeX 才启用 `usetex`，否则用 matplotlib 自带的公式渲染 |
| 6 | `OctupoleTst.py` 一直失败：测试用的是旧版 `vStrt` 格式；Octupole 类没有设置 `Strt2End`/`Rot2LbEnd`，放不进束线；预期值本身前后矛盾 | `02-Tests/OctupoleTst.py`、`BeamLineElement.Octupole.__init__` | 按 Drift 的写法补齐两行；测试改为先建 BeamLine、登记参考粒子，预期值按（线性近似下等价于漂移段的）传输矩阵重新计算 |
| 7 | NumPy 2 下，`plotBeam.py`/`plotextrapolatedBeam.py` 输出 CSV 的 `Twiss` 列写成了 `np.float64(...)`，无法解析 | `01-Code/Beam.py` `setTwiss` | Twiss 值转成 Python float |
| 8 | `03-Scripts`、用户目录里的 9 个脚本和 `02-Tests/decayCHAIN-vis.py` 运行成功时退出码是 1（用法错误时反而是 0），批量运行没法判断成败 | 各脚本最后一行 | 改成 `sys.exit(0)` |
| 9 | 文档小错：README 里的安装指南路径写错；用户目录 README 写的是 `startup.bash`，应为 `user_startup.bash` | `README.md`、`31-UserDirectory/README.md` | 改正 |

## 新增内容

- `02-Tests/runAllTests.bash`：一条命令跑完全部标准测试（按依赖关系排好顺序，约 20 分钟），打印 PASS/FAIL 汇总，日志在 `99-Scratch/TestLogs/`。
- `31-UserDirectory/03-Scripts/scanParameter.py`：参数扫描示例，改参数文件里的一个参数 → 模拟 → 算包络 → 汇总末端束斑和粒子数（CSV + PDF）。
- `.gitignore` 加了 `LaTex.tex`（`LaTexTst.py` 按默认路径在根目录写出的文件）。

## 如何跑全部测试

```bash
source startup.bash                    # 在仓库根目录
bash 02-Tests/runAllTests.bash         # 全部测试
bash 02-Tests/runAllTests.bash DriftTst FocusQuadTst   # 只跑指定的测试
```

## 待办事项

- [x] 写中文代码解读与科研使用指南 `00-Documentation/CodeGuide-CN.md`
- [x] 用 `02-Tests/runAllTests.bash` 做最终验收：61/61 通过（见会话日志）
- [ ] 用户在本地 VS Code 里配置环境，先用**原始**主配置 `LhARABeamLine-Params-LsrDrvn-Gabor.csv` 跑一次基准（`CodeGuide-CN.md` 第 5.1、5.2 节），再跑 `runAllTests.bash`
- [ ] 基准跑通后：复制一份主配置作为自己的参数文件；先和导师确认要用的激光参数（主配置里没写激光参数，用的是代码默认值：10 J、40 fs、焦斑半径 1.5 µm 等）

### 已发现但暂未修改的小问题（需要先确认原作者的本意）

- `03-Scripts/2BDSIM.py`：解析了 `-l <start location>`，但之后没有使用（变量 `iLoc` 未被引用），导出位置固定由 `writeParticleBDSIM(..., 1, True)` 决定。
- 多个脚本的 getopt 长选项定义不完整（如 `"bfile"` 少了 `=`），`--bfile x` 这种写法不起作用；短选项正常。
- `runBeamSim.py` 的 `-z` 用 `bool(arg)` 解析，任何非空字符串（包括 `False`）都当作 True。
- TN 文档（`LhARA-TN-2024-07`）有几处已经过时：说代码“只处理质子”（现在支持 ¹²C⁶⁺、π、μ）；说 `plotextrapolatedBeam.py` 输出 `extrapolatedBeamProgress.pdf`（实际和 `plotBeam.py` 一样写 `BeamProgressionPlot.pdf`）。

- **（2026-10-26 报告中发现，影响结果解读）** 整条束线 91% 的“损失”来自 `BeamLineElement.Transport` 里的纵向截断 `abs(_R[4]) > 5.`（|z| > 5 m），不是物理孔径；全束线传输率和末端能谱主要由它决定。需要和导师确认这个截断的物理含义。
- **（同上）** 偏转段入口的 `LhARA:1:Arc:RPLCswitch:1` 丢掉了 3.5% 的生成粒子（13–17 MeV），原因未查明。
- **（同上）** 线性输运里 x′ = pₓ/p₀，漂移段用参考动量矩阵，Gabor 透镜用按粒子动量计算的物理斜率矩阵，结果是远离 15 MeV 的粒子横向位置被缩放了 p/p₀（比值类的量不受影响）。
- **（同上）** 源模型：能量按全角度能谱抽样，再把方向强制限制在 `rpmax` 内（主配置 5 mrad，只相当于 TNSA 锥内约 0.08% 的质子），所以锥内能谱偏软；能量转换效率定标在 2×10²¹ W/cm² 下已超出适用范围（被截断为 0.5）；质子源尺寸取激光焦斑 r0，而电子鞘层半径是 3.8 µm。

### 科研方面的下一步（详见 `CodeGuide-CN.md` 第 6 节）

1. 复现：跑 LhARA 主配置，和 `41-ReferencePlots/` 对比；比较 `plotBeam` 与 `plotextrapolatedBeam` 的包络。
2. 终端站束流输运优化（原作者 2026-09 的当前重点）：用 `scanParameter.py` 扫末端聚焦四极铁，改进 `EnvelopeOptimisation.py` 的代价函数和优化算法。
3. 根据导师的要求，从指南 6.2 节的方向 B–G 中选题。

## 会话日志

### 2026-10-23 云端会话
- 检查了仓库的可读性：359 个文件全部可读，101 个 `.py` 文件都能编译通过，GitHub 远程访问正常。
- 在最新依赖版本下跑测试，发现并修复了上面表格里的问题 1–4。
- 第一轮全量测试（修复 3、4 之前）：55 个测试脚本里，除了 tag 问题、LaTeX 问题和 OctupoleTst，其余全部通过。
- 第二轮全量测试（修复 1–6 之后）：`02-Tests` 第一阶段 55/55、第二阶段（读取前面生成的数据）5/5 全部通过。
- 第三轮：`03-Scripts` 的 6 个脚本、从 BDSIM 文件读入做束流输运、用户目录（按 README 拷贝后运行）的 4 个脚本，共 12 项全部通过。
- 测试环境：Python 3.13，numpy 2.5.3、pandas 3.0.5、scipy 1.18.1、matplotlib 3.11.2、GitPython 3.2.1、pylorentz；没有 LaTeX。
- 写了中文指南 `00-Documentation/CodeGuide-CN.md`；准直误差（Shift/Tilt 行）和 `scanParameter.py` 都实测有效。
- 最终验收：`bash 02-Tests/runAllTests.bash` 共 61 项（`02-Tests` 里除 `cProfileStats` 外的全部脚本，`decayCHAIN-vis` 带参数运行）全部通过，顺序执行约 18 分钟。
- 之后讨论了本地运行：建议用户在 VS Code 终端里自己跑（Windows 用 WSL），先跑原始主配置作为基准；对话上下文已很长，建议开新对话继续。

### 2026-10-26 云端会话
- 用户在本地跑通前，先请云端跑一次完整模拟并出报告。用主配置 `LhARABeamLine-Params-LsrDrvn-Gabor.csv` 跑了 5 万事件（种子 20261026），报告在 `51-Reports/2026-10-26-Source-GaborLens1/`：`report-CN.md`、`report-EN.md`，分析脚本 `sourceGL1Analysis.py`，图在 `figures/`，结果表在 `results/`。
- 报告重点：激光质子源（TNSA，Te = 10.5 MeV，Kmax = 33.9 MeV）、源处电子的行为与“能否用电子诊断质子源”（结论：部分可行）、第一个 Gabor 透镜（15 MeV 焦距 0.92 m，约 11.3 MeV 出射平行，强色差使投影发射度增长约 480 倍）。
- 新发现的代码问题见上面“已发现但暂未修改的小问题”。
- 数据文件 `99-Scratch/LhARA-LsrDrvn-Gabor.dat`（84 MB）没有提交，用报告第 8 节的命令可以完全重新生成。

