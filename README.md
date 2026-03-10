# POSCARKIT

## Introduction

POSCARKIT 是一个用于处理 VASP POSCAR 文件的工具, 主要用于基于原子占位的建模.
POSCARKIT is a tool for modeling structure POSCAR files, based on Sublattice Occupied Fractions (SOFs).

## Features

- Structure modeling for ordered framework materials.
- Atomic allocation based on Sublattice occupied fractions (SOFs).
- Count Coordinate Numbers (CN) of all pairs of atoms in the supercell and calculate the Nearest Neighbors (NN).
- Slice structure to layers normal to a Miller Index.
- Slice to layers and count CN for each layer.
- Supercell generation along the basis vectors.
- Compare two POSCAR files.
- Merge multiple POSCAR files into one.
- Separate a POSCAR file by different criteria.

## Usage

### Run executable

```shell
# Run it directly
POSCARKIT.exe

# Run it with parameters
POSCARKIT.exe -p <POSCAR.vasp> -c <config.toml>
```

### Run from source

```shell
# Clone the repository
git clone https://github.com/Barabama/POSCARKIT.git
# Build
pip install -e .
# Run
python main.py
```

### Configuration

config.toml 文件用于配置工具行为, 包括超胞生成、原子分配、统计配位、切片等操作的参数.
**可以不进行配置, 对于未配置好的参数程序会请求用户手动输入参数.**

- `name`: 工作名称.
- `poscar`: POSCAR 文件的路径.
- `outdir`: 输出目录.
- `phase`: 指定晶体结构类型（例如 'fcc', 'bcc', 'hcp'）.
- `supercell_factors`: 超胞因子, 沿晶格矢量方向的倍数.
- `slice_direction`: 切片操作的方向, 指定沿晶格矢量方向的切片方向.
- `shuffle_seeds`: 乱序亚晶格位置的随机种子列表.

- `structure_info`: 结构信息参数, 包括结构的单胞结构信息和亚晶格位置的原子分数.
  - Site of Fractions 参数, 用于指定不同亚晶格位置的原子分数.
    （例如, `fcc.1a.sofs`: fcc 结构中亚晶格位置 'a' 的原子分数；
           `bcc.1b.sofs`: bcc 结构中亚晶格位置 'b' 的原子分数.）
  - 单胞结构信息, 这些信息用于生成不同结构类型的单胞 POSCAR 文件.

## CLI 命令行工具

安装后可以使用 `poscarkit` 命令行工具:

### 查看帮助

```bash
poscarkit help
```

### 建模工作流

```bash
poscarkit modeling --poscar POSCAR.vasp --factors 3 3 3 --name my_model
poscarkit modeling --config config.toml --phase fcc --seeds 1 2 3
```

### 配位数统计

```bash
poscarkit countcn --poscar POSCAR.vasp --name my_cn --outdir output/
poscarkit countcn --poscar POSCAR.vasp --cutoff-mult 1.2 --parallel 4
```

### 切片操作

```bash
poscarkit slice --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
```

### 切片并统计配位数

```bash
poscarkit slice-to-countcn --poscar POSCAR.vasp --miller-index 1 1 1 --outdir output/
```

### 超胞生成

```bash
poscarkit supercell --poscar POSCAR.vasp --factors 2 2 2 --outdir output/
```

### 比较 POSCAR 文件

```bash
poscarkit compare --poscar1 POSCAR1.vasp --poscar2 POSCAR2.vasp
```

### 合并 POSCAR 文件

```bash
# 合并两个文件
poscarkit merge --poscars POSCAR1.vasp POSCAR2.vasp --outdir output/

# 合并多个文件
poscarkit merge --poscars POSCAR1.vasp POSCAR2.vasp POSCAR3.vasp --outdir output/
```

### 分离 POSCAR 文件

```bash
# 按 note 分离
poscarkit separate --poscar POSCAR.vasp --key note --outdir output/

# 按 symbol 分离
poscarkit separate --poscar POSCAR.vasp --key symbol --outdir output/
```

## 操作选项 (交互式模式)

运行 `python main.py` 或 `poscarkit_interact` 进入交互式模式:

1. 帮助: 显示帮助信息.
2. 读取配置: 读取 config.toml 文件中的配置.
3. 工作流: 结合超胞生成、分配的工作流.
4. 配位数统计: 统计 POSCAR 文件中每个原子的配位数和最近邻数.
5. 切片: 根据配置对 POSCAR 文件进行切片.
6. 切片并统计配位数: 切片后对每层进行配位数统计.
7. 超胞生成: 根据 POSCAR 文件或配置生成超胞.
8. 比较: 比较两个 POSCAR 文件.
9. 合并: 合并多个 POSCAR 文件.
10. 分离: 根据指定条件分离 POSCAR 文件.

## 项目结构

```
.
├── examples/              # 示例文件
├── old/                   # 旧代码备份
├── src/                   # 核心功能模块
│   ├── cli/               # 命令行接口
│   │   ├── poscarkit.py           # 命令行主程序
│   │   └── poscarkit_interact.py  # 交互式命令行
│   ├── config/            # 配置模块
│   ├── gui/               # GUI 相关资源
│   ├── modeling/          # 建模核心模块
│   │   ├── __init__.py
│   │   ├── base.py        # 基础结构类
│   │   ├── countcn.py     # 配位数统计
│   │   ├── model.py       # 模型相关
│   │   ├── slice.py       # 切片操作
│   │   └── supercell.py   # 超胞生成
│   ├── utils/             # 通用工具
│   └── workflow/          # 工作流模块
│       ├── modeling.py        # 建模工作流
│       └── slice_to_countcn.py # 切片到配位数统计工作流
├── tests/                 # 测试代码
│   ├── cli/               # CLI 测试
│   ├── modeling/          # 建模模块测试
│   ├── utils/             # 工具测试
│   └── workflow/          # 工作流测试
├── main.py                # 主程序入口
├── pyproject.toml         # 项目配置
└── README.md              # 项目文档
```

## 工作流

### 建模工作流 (Modeling Workflow)

建模工作流包含两个主要步骤: 
1. 超胞生成: 根据指定的超胞因子生成超胞
2. 原子分配: 根据结构信息和随机种子分配原子

该工作流实现在 `src/workflow/modeling.py` 模块中.

### 切片到配位数统计工作流 (Slice to Count CN Workflow)

此工作流结合了切片和配位数统计功能: 
1. 按指定晶面指数切片
2. 对每个切片层分别进行配位数统计

该工作流实现在 `src/workflow/slice_to_countcn.py` 模块中.

## 核心功能

### 1. 超胞生成 (Supercell Generation)

根据输入的 POSCAR 文件和超胞因子, 生成指定大小的超胞结构.

### 2. 原子分配 (Atomic Allocation)

根据配置文件中的亚晶格原子分数, 对超胞中的原子进行随机分配, 生成具有指定成分的结构.

### 3. 配位数统计 (Coordination Number Counting)

统计结构中每个原子的配位数和最近邻原子对的数量, 生成配位数分布直方图.

### 4. 切片操作 (Slicing)

将结构沿指定方向切片, 生成平行于指定晶面的原子层结构, 并可选择对每层进行配位数统计.

### 5. POSCAR 合并与分离 (Merge & Separate)

- 合并多个 POSCAR 文件为一个
- 根据指定条件分离一个 POSCAR 文件为多个

## 测试

本项目包含完整的单元测试, 可以通过以下方式运行: 

### 运行所有测试

```bash
python -m unittest discover tests
```

### 运行特定模块的测试

```bash
python -m unittest tests.modeling.test_simple_poscar
python -m unittest tests.modeling.test_supercell
python -m unittest tests.workflow.test_modeling
python -m unittest tests.workflow.test_slice_to_countcn
```

## 编译

- 用Nuitka打包

```Shell
# 安装Nuitka
cd poscarkit
pip install nuitka

# nuitka --standalone --onefile --output-dir=dist --jobs4 --lto=yes `
# --follow-imports --enable-plugin=no-qt --include-package=pandas `
# --nofollow-import-to=matplotlib.tests --nofollow-import-to=pandas.tests `
# --nofollow-import-to=pytest --nofollow-import-to=setuptools.tests `
# --windows-icon-from-ico="src/gui/poscarkit.ico" `
# --output-filename=poscarkit-0.9.2.exe `
# --file-version=0.9.2 `
# --copyright="(C) 2025 MCMF, Fuzhou University" `
# main.py
```
