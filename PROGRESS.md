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

## 待办事项

- [ ] 没装 LaTeX 时 `decayCHAIN.py`、`ExpSourceTst.py`、`ExpSourceTst_rms_KE.py` 作图失败（`usetex=True`）→ 改成有 LaTeX 才启用
- [ ] `OctupoleTst.py` 失败（原作者笔记里一直记为 Fails）→ 待查
- [ ] 全部修完后重跑 `02-Tests` 全量测试
- [ ] 测试 `03-Scripts/` 和 `31-UserDirectory/` 里的脚本（需要命令行参数）
- [ ] 写中文代码解读与科研使用指南 `00-Documentation/CodeGuide-CN.md`

## 会话日志

### 2026-10-23 云端会话
- 检查了仓库的可读性：359 个文件全部可读，101 个 `.py` 文件都能编译通过，GitHub 远程访问正常。
- 在最新依赖版本下跑测试，发现并修复了上面表格里的问题 1–4。
- 第一轮全量测试（修复 3、4 之前）：55 个测试脚本里，除了 tag 问题、LaTeX 问题和 OctupoleTst，其余全部通过。
