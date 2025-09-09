# POSCARKIT

## 介绍

POSCARKIT 是一个用于处理 VASP POSCAR 文件的工具，主要用于基于原子占位的建模。

## 使用说明

以下是对于可执行文件的使用说明：

### 1. 打开方式

- 支持双击打开。
- 支持POSCAR拖拽使用该文件打开。
- 支持命令行参数打开。
`POSCARKIT.exe <POSCAR.vasp>`
`POSCARKIT.exe -f <POSCAR.vasp> -c <choice>`

P.S. 不指定参数则默认读取 config.toml 文件的配置。

### 2. 操作选项

1. 帮助：显示帮助信息。
2. 读取配置：读取 config.toml 文件中的配置。
3. 工作流：结合超胞生成、分配的工作流。
4. 超胞生成：根据 POSCAR 文件或配置生成超胞。
5. 分配：根据配置对 POSCAR 文件中的原子进行随机洗牌和分配。
6. 统计配位：统计 POSCAR 文件中每个原子的配位数和最近邻数。
7. 切片：根据配置对 POSCAR 文件进行切片。

### 3. config.toml 配置

config.toml 文件用于配置工具行为，包括超胞生成、原子分配、统计配位、切片等操作的参数。
**可以不进行配置，对于未配置好的参数程序会请求用户手动输入参数。**

- `Filepath`：POSCAR 文件的路径。
- `Outdir`：输出文件夹的路径。
- `SupercellFactors`：超胞生成的因子，用于指定沿晶格矢量方向的倍数。
- `Structure`：指定结构类型（例如 "fcc", "bcc", "hcp"）。
- `ShuffleSeeds`：洗牌操作使用的随机种子列表。
- `SliceDirection`：切片操作的方向，指定沿晶格矢量方向的切片方向。

- `StructureInfo`：结构信息参数，包括结构的单胞结构信息和亚晶格位置的原子分数。
  - Site of Fractions 参数，用于指定不同亚晶格位置的原子分数。
    （例如, `FCC.1a.sofs`：FCC 结构中亚晶格位置 'a' 的原子分数；
           `BCC.1b.sofs`：BCC 结构中亚晶格位置 'b' 的原子分数。）
  - 单胞结构信息，这些信息用于生成不同结构类型的单胞 POSCAR 文件。

## 安装

```bash
git clone https://gitee.com/wubo-movers/poscarkit.git  # 克隆仓库
cd poscarkit
# 选择一个python虚拟环境
pip install -r requirements.txt  # 安装依赖
python poscarkit.py  # 运行程序
```

## 功能特性

- 超胞生成：根据 POSCAR 文件或配置生成超胞
- 原子分配：对原子进行随机洗牌和分配
- 统计配位：统计每个原子的配位数和最近邻数
- 切片：根据配置对 POSCAR 文件进行切片
- 工作流：包含建模工作流和切片与配位数统计工作流

## 项目结构

```
.
├── PoscarTools/           # 核心功能模块
│   ├── AtomAllocate.py    # 原子分配模块
│   ├── AtomCountCN.py     # 统计配位数模块
│   ├── AtomSlice.py       # 切片操作模块
│   ├── AtomSupercell.py   # 超胞生成模块
│   ├── SimplePoscar.py    # POSCAR 文件基础解析模块
│   └── Utils.py           # 通用工具函数
├── workflows/             # 工作流模块
│   ├── modeling.py        # 建模工作流（超胞生成+原子分配）
│   └── sliceandcountcn.py # 切片和配位数统计工作流
├── tests/                 # 测试代码
│   ├── PoscarTools/       # 各核心模块的单元测试
│   └── workflows/         # 工作流模块的单元测试
├── examples/              # 示例文件
├── poscarkit.py           # 主程序入口
└── config.toml            # 配置文件模板
```

## 工作流

### 建模工作流 (Modeling Workflow)

建模工作流包含两个主要步骤：
1. 超胞生成：根据指定的超胞因子生成超胞
2. 原子分配：根据结构信息和随机种子分配原子

该工作流已从主程序移到 `workflows/modeling.py` 模块中，通过 `run_modeling_workflow` 函数实现。

### 切片和配位数统计工作流 (Slice and Count CN Workflow)

此工作流结合了切片和配位数统计功能：
1. 按指定晶面指数切片
2. 对每个切片层分别进行配位数统计

该工作流实现在 `workflows/sliceandcountcn.py` 模块中。

## 测试

本项目包含完整的单元测试，可以通过以下方式运行：

### 运行所有测试

```bash
python -m unittest tests.test_suite
```

或者：

```bash
python tests/test_suite.py
```

### 运行特定模块的测试

```bash
python -m unittest tests.PoscarTools.test_atom_supercell
python -m unittest tests.PoscarTools.test_atom_allocate
python -m unittest tests.PoscarTools.test_atom_slice
python -m unittest tests.workflows.test_modeling
python -m unittest tests.workflows.test_sliceandcountcn
```

### 测试覆盖率

目前测试覆盖了以下核心功能：

1. 超胞生成（AtomSupercell）
2. 原子分配（AtomAllocate）
3. 切片操作（AtomSlice）
4. 工作流功能（workflows）

## 编译

- 用Nuitka和UPX打包

```Shell
# 安装UPX
curl -O https://github.com/upx/upx/releases/download/v5.0.0-win64/upx-5.0.0-win64.zip
unzip upx-5.0.0-win64.zip
cd upx-5.0.0-win64
setx PATH "%PATH%;$PWD.Path"
upx --version

# 安装Nuitka
cd poscarkit
pip install nuitka
nuitka poscarkit.py --standalone --onefile --output-dir=dist --remove-output `
--windows-icon-from-ico="icon.ico" `
--enable-plugin=upx --upx-binary="D:\Programs\upx-5.0.2-win64\upx.exe" `
--enable-plugin=tk-inter `
--follow-imports `
```
