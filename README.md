# POSCARKIT

**A toolkit for modeling VASP POSCAR files based on Sublattice Occupying Fractions (SOFs).**
**基于亚晶格占位分数 (SOFs) 的 VASP POSCAR 结构建模工具包.**

[![Version](https://img.shields.io/badge/version-0.10.6-blue)](./pyproject.toml)
[![Python](https://img.shields.io/badge/python-≥3.10-blue)](./pyproject.toml)
[![License](https://img.shields.io/badge/license-MIT-green)](./LICENSE)

---

## Features · 功能

| Feature · 功能 | Description · 描述 |
|---|---|
| Modeling · 建模 | Generate supercell and allocate atoms based on SOFs · 基于 SOFs 生成超胞并分配原子 |
| Import to Model · 导入建模 | Import SOFs from CSV/XLSX and run modeling · 从 CSV/XLSX 导入 SOFs 并运行建模 |
| Count CN · 配位数统计 | Coordination number & nearest-neighbor pair counting · 配位数和最近邻原子对统计 |
| Slice · 切片 | Slice structure normal to a Miller index · 沿晶面指数切片 |
| Slice to CountCN · 切片配位 | Slice then count CN per layer · 切片后逐层统计配位数 |
| Surface · 表面 | Generate asymmetric surface slabs from bulk POSCAR with vacuum and selective dynamics · 从体相 POSCAR 生成非对称表面 slab 模型（含真空层和选择性动力学约束） |
| Thermo · 热力学 | Calculate configurational entropy (Sconf) and Gibbs free energy (DeltaG) from SOF data + TDB · 从 SOF 数据和 TDB 热力学数据库计算构型熵和 Gibbs 自由能 |
| Supercell · 超胞 | Expand unit cell along basis vectors · 沿基矢方向扩展超胞 |
| Compare · 比较 | Diff two POSCAR files · 比较两个 POSCAR 文件 |
| Merge · 合并 | Combine multiple POSCARs into one · 合并多个 POSCAR 文件 |
| Separate · 分离 | Split a POSCAR by element, coordinate, or group · 按元素/坐标/分组分离 POSCAR |

---

## Quick Start · 快速开始

### Run executable · 运行可执行文件

```bash
POSCARKIT.exe                           # launch GUI · 启动图形界面
POSCARKIT.exe -p POSCAR.vasp -c config.toml
```

### Run from source · 从源码运行

```bash
git clone https://github.com/Barabama/POSCARKIT.git
cd poscarkit
pip install -e .
python main.py                          # no args → GUI · 无参数→图形界面
python main.py modeling --help          # with args → CLI · 带参数→命令行
```

---

## GUI · 图形界面

`python main.py` (no arguments) opens a Tkinter GUI window with sidebar navigation.
`python main.py`（无参数）打开 Tkinter 图形界面，左侧导航栏切换功能。

- **Sidebar · 侧边栏**: Modeling / Import to Model / Count CN / Slice / Slice to CountCN / Surface / Thermo / Supercell / Compare / Merge / Separate
- **Form area · 表单区**: each function shows a parameter form pre-filled from `config.toml`
- **SOF editor · 占位编辑器**: Modeling form includes an inline SOF table (add/remove elements, real-time sum check)
- **Log area · 日志区**: real-time logging output at the bottom
- **Run in background · 后台运行**: tasks execute in a separate thread, UI stays responsive
- **Save Config · 保存配置**: button on the Config form to persist all settings to `config.toml`
- **Load from config · 从配置加载**: button on each form to reload the latest config values into the form
- Scrollable form area with draggable form/log divider · 表单区可滚动，表单/日志分隔线可拖拽

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

### Thermo · 热力学

```bash
poscarkit thermo --data sof_data.xlsx --tdb database.TDB --name my_thermo
poscarkit thermo --data sof_data.xlsx --tdb database.TDB --name my_thermo --outdir output/
```

### Import to Model · 导入建模

```bash
poscarkit import-to-model --csv tc_exps.csv --phase fcc -t 473 873 1273
poscarkit import-to-model --csv pandat-data.csv --phase fcc -t 873 -o print
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

### Surface · 表面

```bash
poscarkit surface POSCAR.vasp --miller 1 1 1 --layers 3 --vacuum 15.0 --name my_slab
poscarkit surface POSCAR.vasp --miller 1 1 0 --layers 4 --fix-layers 2 --fix-z-only
```

Outputs N and N+1 layer slabs per termination with Selective Dynamics (F F F or T T F),
vacuum (bottom 2 Å + top), and a summary CSV with composition deviation, dipole moment, etc.
输出每个终止面的 N 和 N+1 层 slab，包含真空层、选择性动力学约束和汇总 CSV。

## Interactive Mode · 交互模式

Launch the interactive numbered menu:

```bash
python -m src.cli.poscarkit_interact
```

```
 1) Help                5) Slice to layers    9)  Merge
 2) Read config         6) Slice to CountCN  10) Separate
 3) Run Modeling        7) Make Supercell    11) Import to Model
 4) Count CN            8) Compare           12) Thermo Analysis
                                          13) Surface Slab
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
│   ├── io/                      # SOF data import · SOF 数据导入
│   │   ├── ir.py                #    intermediate representation · 中间表示
│   │   └── readers/             #    format-specific readers · 格式读取器
│   ├── modeling/                # core modeling · 核心建模
│   │   ├── base.py              #    Atom, Struct, SimplePoscar
│   │   ├── countcn.py           #    CN counter with KDTree & ASE backends
│   │   ├── model.py             #    ModelStruct — shuffle & SQS allocation
│   │   ├── slice.py             #    Slicer — Miller-index layer slicing
│   │   ├── surface.py           #    SurfaceBuilder — slab generation · 表面 slab 建模
│   │   └── supercell.py         #    supercell generation · 超胞生成
│   ├── thermo/                  # thermodynamics · 热力学计算
│   │   ├── tdb.py               #    TDB parser · TDB 解析器
│   │   ├── sconf.py             #    configurational entropy · 构型熵
│   │   ├── deltaG.py            #    Gibbs free energy · Gibbs 自由能
│   │   └── plot.py              #    thermo plots · 热力学绘图
│   ├── utils/                   # utilities · 工具
│   └── workflow/                # high-level workflows · 高层工作流
│       ├── import_to_model.py   #    import SOFs from CSV/XLSX and run modeling
│       ├── modeling.py          #    supercell + allocation · 超胞+分配
│       ├── slice_to_countcn.py  #    slice + CN per layer · 切片+逐层配位
│       └── thermo.py            #    Sconf + DeltaG pipeline · 热力学管道
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
python -m unittest tests.modeling.test_slice
python -m unittest tests.modeling.test_supercell
python -m unittest tests.modeling.test_surface
python -m unittest tests.workflow.test_modeling
```

---

## Build · 编译

Packaged with Nuitka · 用 Nuitka 打包:

```bash
pip install nuitka
nuitka --standalone --onefile --output-dir=dist --jobs=4 --lto=yes \
    --follow-imports --enable-plugins=tk-inter --include-package=pandas  \
    --nofollow-import-to=matplotlib.tests --nofollow-import-to=pandas.tests \
    --nofollow-import-to=pytest --nofollow-import-to=setuptools.tests \
    --windows-icon-from-ico="src/gui/poscarkit.ico" --windows-console-mode=disable \
    main.py
```

---

## License · 许可证

MIT © 2025 MCMF, Fuzhou University · MCMF
