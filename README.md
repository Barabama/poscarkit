# POSCARKIT

**A toolkit for modeling VASP POSCAR files based on Sublattice Occupying Fractions (SOFs).**
**基于亚晶格占位分数 (SOFs) 的 VASP POSCAR 结构建模工具包.**

[![Version](https://img.shields.io/badge/version-0.10.1-blue)](./pyproject.toml)
[![Python](https://img.shields.io/badge/python-≥3.10-blue)](./pyproject.toml)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)

---

## Features · 功能

| Feature · 功能 | Description · 描述 |
|---|---|
| Modeling · 建模 | Generate supercell and allocate atoms based on SOFs · 基于 SOFs 生成超胞并分配原子 |
| Count CN · 配位数统计 | Coordination number & nearest-neighbor pair counting · 配位数和最近邻原子对统计 |
| PBC support · 周期性 | Optional periodic boundary conditions for CN counting · 可选的周期性边界条件配位数统计 |
| Slice · 切片 | Slice structure normal to a Miller index · 沿晶面指数切片 |
| Slice + CN · 切片配位 | Slice then count CN per layer · 切片后逐层统计配位数 |
| Supercell · 超胞 | Expand unit cell along basis vectors · 沿基矢方向扩展超胞 |
| Compare · 比较 | Diff two POSCAR files · 比较两个 POSCAR 文件 |
| Merge · 合并 | Combine multiple POSCARs into one · 合并多个 POSCAR 文件 |
| Separate · 分离 | Split a POSCAR by element, coordinate, or group · 按元素/坐标/分组分离 POSCAR |

---

## Quick Start · 快速开始

### Run executable · 运行可执行文件

```bash
POSCARKIT.exe                           # interactive mode · 交互模式
POSCARKIT.exe -p POSCAR.vasp -c config.toml
```

### Run from source · 从源码运行

```bash
git clone https://github.com/Barabama/POSCARKIT.git
cd poscarkit
pip install -e .
python main.py                          # no args → interactive · 无参数→交互模式
```

---

## CLI Commands · 命令行

```
poscarkit help                          # show help · 显示帮助
```

### Modeling · 建模

```bash
poscarkit modeling --poscar POSCAR.vasp --factors 3 3 3 --name my_model
poscarkit modeling --config config.toml --phase fcc --seeds 1 2 3
poscarkit modeling --config config.toml --phase fcc --enable-sqs --batch-size 4
```

### Count CN · 配位数统计

```bash
poscarkit countcn --poscar POSCAR.vasp --name my_cn
poscarkit countcn --poscar POSCAR.vasp --cutoff-mult 1.2 --parallel 4
poscarkit countcn --poscar POSCAR.vasp --by-ase              # use ASE backend · 使用 ASE 后端
poscarkit countcn --poscar unitcell.vasp --pbc               # with PBC · 启用周期性边界
```

### Slice · 切片

```bash
poscarkit slice --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
```

### Slice to Count CN · 切片后统计配位

```bash
poscarkit slice-to-countcn --poscar POSCAR.vasp --miller-index 1 1 1
poscarkit slice-to-countcn --poscar POSCAR.vasp --miller-index 1 1 1 --pbc
```

### Supercell · 超胞生成

```bash
poscarkit supercell --poscar POSCAR.vasp --factors 2 2 2 --outdir output/
poscarkit supercell --poscar POSCAR.vasp --factors 2 2 2 --by-ase
```

### Compare · 比较

```bash
poscarkit compare --poscar1 POSCAR1.vasp --poscar2 POSCAR2.vasp
```

### Merge · 合并

```bash
poscarkit merge --poscars POSCAR1.vasp POSCAR2.vasp --outdir output/
poscarkit merge --poscars A.vasp B.vasp C.vasp --outdir output/
```

### Separate · 分离

```bash
poscarkit separate --poscar POSCAR.vasp --key note --outdir output/
poscarkit separate --poscar POSCAR.vasp --key symbol --outdir output/
```

## Interactive Mode · 交互模式

`python main.py` (no arguments) enters a numbered menu · `python main.py`（无参数）进入菜单:

```
 1) Help                5) Slice to layers    8)  Compare
 2) Read config         6) Slice to CountCN   9)  Merge
 3) Run Modeling        7) Make Supercell    10)  Separate
 4) Count CN
```

---

## Configuration · 配置

The `config.toml` file defines phases, sublattice geometries, and SOF parameters.
`config.toml` 文件定义晶体相、亚晶格几何和 SOF 参数。

Key parameters · 关键参数:

| Parameter · 参数 | Description · 描述 |
|---|---|
| `name` | Work name · 工作名称 |
| `poscar` | Input POSCAR path · 输入文件路径 |
| `outdir` | Output directory · 输出目录 |
| `phase` | Crystal structure: `fcc`, `bcc`, `hcp` · 晶体结构类型 |
| `supercell_factors` | Expansion factors e.g. `[3, 3, 3]` · 超胞因子 |
| `shuffle_seeds` | Random seeds for reproducibility · 随机种子 |
| `batch_size` | Number of parallel modeling batches · 并行批次数 |
| `enable_sqs` | Enable SQS generation · 启用 SQS 生成 |
| `cutoff_mult` | CN cutoff multiplier · 配位截断倍数 |
| `by_ase` | Use ASE for CN / supercell · 使用 ASE 后端 |
| `pbc` | Periodic boundary conditions for CN · 配位统计周期边界 |

SOFs are defined per phase per sublattice · SOFs 按相/亚晶格定义:

```toml
[fcc.1a.sofs]
V = 1.0

[fcc.3c.sofs]
Co = 0.4444
Ni = 0.4444
V = 0.1112
```

Custom phases can be added · 可添加自定义相 (see `config.toml` Advanced Usage).

---

## Project Structure · 项目结构

```
.
├── examples/                    # sample POSCARs · 示例文件
├── old/                         # legacy scripts · 旧版脚本
├── src/
│   ├── cli/                     # command-line interface · 命令行接口
│   │   ├── poscarkit.py         #    CLI entry · 命令行入口
│   │   └── poscarkit_interact.py#    interactive menu · 交互菜单
│   ├── config/                  # metadata & default config · 元数据与默认配置
│   ├── modeling/                # core modeling · 核心建模
│   │   ├── base.py              #    Atom, Struct, SimplePoscar
│   │   ├── countcn.py           #    CN counter with KDTree & ASE backends
│   │   ├── model.py             #    ModelStruct — shuffle & SQS allocation
│   │   ├── slice.py             #    Slicer — Miller-index layer slicing
│   │   └── supercell.py         #    supercell generation · 超胞生成
│   ├── utils/                   # utilities · 工具
│   └── workflow/                # high-level workflows · 高层工作流
│       ├── modeling.py          #    supercell + allocation · 超胞+分配
│       └── slice_to_countcn.py  #    slice + CN per layer · 切片+逐层配位
├── tests/                       # unit tests · 单元测试
├── main.py                      # program entry · 程序入口
├── pyproject.toml               # project metadata & deps · 项目元数据
├── config.toml                  # user configuration · 用户配置
└── README.md
```

---

## Testing · 测试

```bash
python -m unittest discover tests             # all tests · 全部测试
python -m unittest tests.modeling.test_simple_poscar
python -m unittest tests.modeling.test_countcn
python -m unittest tests.workflow.test_modeling
```

---

## Build · 编译

Packaged with Nuitka · 用 Nuitka 打包:

```bash
pip install nuitka
nuitka --standalone --onefile --output-dir=dist --jobs=4 --lto=yes \
    --follow-imports --enable-plugin=no-qt --include-package=pandas \
    --nofollow-import-to=matplotlib.tests --nofollow-import-to=pandas.tests \
    --nofollow-import-to=pytest --nofollow-import-to=setuptools.tests \
    --windows-icon-from-ico="src/gui/poscarkit.ico" \
    --output-filename=poscarkit-0.10.1.exe \
    --file-version=0.10.1 \
    --copyright="(C) 2025 MCMF, Fuzhou University" \
    main.py
```

---

## License · 许可证

MIT © 2025 MCMF, Fuzhou University · MCMF
