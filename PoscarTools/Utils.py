# utils.py

color_map = {"Li": "pink",
             "B": "sienna",
             "Mg": "lawngreen",
             "Al": "silver",
             "Ti": "gray",
             "V": "yellow",
             "Cr": "red",
             "Mn": "purple",
             "Fe": "black",
             "Co": "blue",
             "Ni": "green",
             "Cu": "peru",
             "Zn": "deepskyblue",
             "Y": "olive",
             "Zr": "orangered",
             "Nb": "darkviolet",
             "Mo": "orange",
             "Ag": "silver",
             "Sn": "slategray",
             "Hf": "deepskyblue",
             "Ta": "darkgoldenrod",
             "W": "gold",
             "Re": "mediumorchid",
             }

default_config = """# config.toml

# P.S. No need to configure everything. Just write down SOFs data and
#      run the program and it will ask you to input the information.
# P.S. 无需配置所有内容。只需写下占位分数并运行程序，它就会要求您输入信息。

# ====== Operations on POSCAR, commenting out what won't be used ======
# ====== 程序操作使用以下配置参数, 不使用的参数请注释掉以免被程序读取 ======
# Filepath of POSCAR  # POSCAR 文件的路径
# Filepath = "./examples/CoNiV.vasp"

# Output directory  # 输出目录
# Outdir = "./out/"

# Used for 'Supercell', 'Allocation'  # 指定超胞因子, 沿晶格矢量方向的倍数
# SupercellFactors = [30, 30, 30]

# Used for 'Allocation'  # 指定晶体结构类型, 不区分大小写
# Structure = "fcc"

# Used for 'Allocation'  # 乱序亚晶格位置的随机种子列表
# ShuffleSeeds = [7, 42, 83]

# Used for 'Slice'  # 指定切片操作的方向，基于晶格矢量方向, 如 [001] z 轴
# SliceDirection = [0, 0, 1]  # Direction of slicing Used for 'Slice'

# ====== Site of Fractions parameters, used for 'Allocation' ======
# ====== 晶体结构对应亚晶格的占位分数的配置参数, 用于 '分配' 操作 ======
# Sum of sofs must be close to 1 # 分数的和必须接近于1
[FCC.1a.sofs]
V = 1.0

[FCC.3c.sofs]
Co = 4.444E-1
Ni = 0.4444
V = 0.1112

[BCC.1a.sofs]

[BCC.1b.sofs]

[HCP.2a.sofs]

[HCP.6c.sofs]


# ====== ! Unit cell structure parameters, be cautious when modifying ! ======
# ====== ! 原型单胞的晶体结构参数, 从这里生成对应结构的单胞, 修改时请小心 ! ======
# It is wise to give the Volume relaxed lattice parameter
# when finished the R7 step (ISIF=7, NSW=10) in INCAR to run VASP
# when the HEA-POSCAR established based on SOFs and small supercell.
# Then we can give considerable reliable ruler tag.
[FCC]
cell = [3.774, 3.774, 3.774]
1a.atoms = ["Au", [[0, 0, 0]]]
3c.atoms = ["Cu", [[0, 0.5, 0.5], [0.5, 0, 0.5], [0.5, 0.5, 0]]]

[BCC]
cell = [2.7, 2.7, 2.7]
1a.atoms = ["Al", [[0, 0, 0]]]
1b.atoms = ["Ni", [[0.5, 0.5, 0.5]]]

[HCP]
cell = [[2.6525, -4.593, 0], [2.6525, 4.593, 0], [0, 0, 4.242]]
2a.atoms = ["Sn", [[0.3333, 0.6667, 0.25], [0.6667, 0.3333, 0.75]]]
6c.atoms = ["Ni",[[0.1596, 0.8404, 0.75], [0.1596, 0.3192, 0.75],
                  [0.3192, 0.1596, 0.25], [0.6808, 0.8404, 0.75],
                  [0.8404, 0.6808, 0.25], [0.8404, 0.1596, 0.75]]]
"""
