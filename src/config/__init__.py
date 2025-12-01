# src/config/__init__.py

# Application metadata
DEVELOPER = "Author: Gao Min-Liang, Qiao Yang, Yang Su-wen, Wu Bo*, et al."
VERSION = "Version: 0.9.0 | © 2025 MCMF, Fuzhou University"
CONTACT = "Email: wubo@fzu.edu.cn | Phone: +86 130 2381 9517"

DEFAULT_CONFIG = """# config.toml

# P.S. No need to configure everything. Just write down SOFs data and
#      run the program and it will ask you to input the information.
# P.S. 无需配置所有内容, 只需写下占位分数并运行程序，它就会要求您输入信息

# ====== Operations on POSCAR, commenting out what won't be used ======
# ====== 程序操作使用以下配置参数, 不使用的参数请注释掉以免被程序读取 ======

# Work name         # 工作名称
# name = "MyWork"

# POSCAR file path  # POSCAR 文件的路径
# poscar = "./examples/CoNiV.vasp"

# Output directory  # 输出目录
# outdir = "./output/"

# Type of phase     # 指定晶体结构类型
# phase = "fcc"

# Supercell factors along lattice vector, default is [3, 3, 3]
# 超胞因子, 沿晶格矢量方向的倍数, 默认为[3, 3, 3]
# supercell_factors = [3, 3, 3]

# Slice direction along lattice vector, default is [0, 0, 1] z-axis
# 切片方向，基于晶格矢量方向, 默认为[0, 0, 1]即z轴
# slice_direction = [0, 0, 1]

# Shuffle seeds for modeling, for reproduction, default not used
# 乱序亚晶格位置的随机种子列表, 用于复现, 默认不使用
# shuffle_seeds = [7, 42, 83]

# Batch size for modeling   # 建模的批次大小
# batch_size = 3

# Enable SQS        # 启用SQS
# enable_sqs = false



# ====== Sublattices Occupied Fractions parameters ======
# ====== 各个晶体结构所对应亚晶格的占位分数的配置参数 ======
# If the sum of sofs is not close to 1, it will ask for normalization.
# 如果占位分数和总和不接近1, 它会询问归一化.

[bcc.1a.sofs]

[bcc.1b.sofs]

[fcc.1a.sofs]
V = 1.0

[fcc.3c.sofs]
Co = 4.444E-1
Ni = 0.4444
V = 0.1112

[hcp.2a.sofs]

[hcp.6c.sofs]


# ====== ! Unit cell structure parameters, be cautious when modifying ! ======
# ====== ! 原型单胞的晶体结构参数, 从这里生成对应结构的单胞, 修改时请小心 ! ======
# It is wise to give the Volume relaxed lattice parameter
# when finished the R7 step (ISIF=7, NSW=10) in INCAR to run VASP
# when the HEA-POSCAR established based on SOFs and small supercell.
# Then we can give considerable reliable ruler tag.
[bcc]
cell = [2.7, 2.7, 2.7]
1a.atoms = ["Al", [[0, 0, 0]]]
1b.atoms = ["Ni", [[0.5, 0.5, 0.5]]]

[fcc]
cell = [3.774, 3.774, 3.774]
1a.atoms = ["Au", [[0, 0, 0]]]
3c.atoms = ["Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]]

[hcp]
cell = [[2.6525, -4.593, 0], [2.6525, 4.593, 0], [0, 0, 4.242]]
2a.atoms = ["Sn", [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]]]
6c.atoms = ["Ni", [[0.1596, 0.8404, 0.75], [0.1596, 0.3192, 0.75],
                  [0.3192, 0.1596, 0.25], [0.6808, 0.8404, 0.75],
                  [0.8404, 0.6808, 0.25], [0.8404, 0.1596, 0.75]]]


# ====== Advanced Usage ======
# It is allowed to add custom structures

# phase = "MyPhase"
# [MyPhase]
# cell = [[2, 2, 0],[2, 2, 0],[0, 0, 3.8]]
# subl1.atoms=["Au", [[0, 0, 0]]]
# subl2.atoms=["Cu", [[0.5, 0.5, 0.5]]]
# [MyPhase.subl1.sofs]
# Au = 1.0
# [MyPhase.subl2.sofs]
# Cu = 1.0
"""
