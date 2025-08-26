# POSCARKIT

## 介绍

POSCARKIT 是一个用于处理 VASP POSCAR 文件的工具，支持多种操作，主要用于基于原子占位的
POSCAR 建模，还支持统计配位、切片等功能。

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

```shell
git clone https://gitee.com/wubo-movers/poscarkit.git  # 克隆仓库
cd poscarkit
# 选择一个python虚拟环境
python -m pip install -r requirements.txt  # 安装依赖
python -m poscarkit.py  # 运行程序
```

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
