# POSCARKIT

## 介绍

POSCARKIT 是一个用于处理 VASP POSCAR 文件的工具，支持多种操作，
包括超胞生成、原子洗牌、原子分配、切片等。

## 使用说明

以下是对于可执行文件的使用说明：

### 1. 打开方式

- 支持双击打开。
- 支持POSCAR拖拽使用该文件打开。
- 支持命令行参数打开。
`POSCARKIT.exe <POSCAR.vasp>`
`POSCARKIT.exe -f <POSCAR.vasp> -c <choice>`

PS. 不指定参数则默认读取 config.toml 文件的配置。

### 2. 操作选项

1. 读取配置：读取 config.toml 文件中的配置。
2. 超胞生成：根据配置生成超胞文件。
3. 切片：根据配置对 POSCAR 文件进行切片。
4. 洗牌：根据配置对 POSCAR 文件中的原子进行洗牌。
5. 分配：根据配置对 POSCAR 文件中的原子进行分配。
6. 工作流：结合超胞生成、洗牌和分配操作。
7. 统计配位数：统计 POSCAR 文件中每个原子的配位数。
8. Exit：Ctrl+C 退出工具。

### 3. config.toml 配置

config.toml 文件用于配置工具行为，包括超胞生成、切片、洗牌和分配等操作的参数。
**可以不进行配置，在程序中手动输入参数。**

- `FilePath`：POSCAR 文件的路径。
- `SupercellFactors`：超胞生成的因子，用于指定沿晶格矢量方向的倍数。
- `Structure`：指定结构类型（例如 "fcc", "bcc", "hcp"）。
- `ShuffleSeeds`：洗牌操作使用的随机种子列表。
- `Shuffle`：是否在分配前进行洗牌操作。
- `SliceDirection`：切片操作的方向，指定沿晶格矢量方向的切片方向。

- Site of Fractions 参数，这些参数用于指定不同 Wyckoff 位置的原子分数
  （例如，`FCC.1a.sofs`：FCC 结构中 Wyckoff 位置 'a' 的原子分数；
  `BCC.1b.sofs`：BCC 结构中 Wyckoff 位置 'b' 的原子分数。）。

- 单元胞结构参数
这些参数用于指定不同结构类型的单元胞参数。（晶体常数基本不使用）
