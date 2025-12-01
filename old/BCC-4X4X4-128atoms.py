from typing import Match
from numpy.core.records import array
import pandas as pd
import numpy as np
from vaspy.atomco import PosCar
import os
import math
from math import floor, modf


############################################ 基础配置信息 #####################################################
workPath = 'E:\Desktop\FCC-BCC-HCP-POSCAR-TemplatePython\BCC\BCC-4X4X4-128atoms'

############################################## 元素个数 #######################################################
AtomNumber = 5

########################################## 1a 亚晶格下的占位分数 ################################################
AtomicOccupationScoreDic = {
    'Al':0.1974,
    'Co':0.3908,
    'Cr':0.2423,
    'Fe':0.0435,
    'Ni':0.1259
}

############################################# POSCAR 基矢系数 ##################################################
poscarBases = np.array([
    [11.39999962,0.00000000,0.00000000],
    [0.00000000,11.39999962,0.00000000],
    [0.00000000,0.00000000,11.39999962]
])

## 亚晶格类型
sublattice1 = '1a'
sublattice2 = '1b'
sublattice3 = '1a+1b'

## 4X4X4 超胞原子总数
atomSum = 128

## 每类亚晶格的原子数
AtomNumberof1 = 64
AtomNumberof2 = atomSum - AtomNumberof1

## 超胞数量类型
strOcc = '4X4X4'

## 某一原子的原子总数
singleAtom = AtomNumberof1

AtomNumberofSingleAtom = math.ceil(atomSum/AtomNumber)

################################################# 计算过程 ######################################################

## 1a亚晶格各个原子的个数
AtomNumberof1Dic = {}
AtomNumberof1List = []

## 1b亚晶格各个原子的个数
AtomNumberof2Dic = {}
AtomNumberof2list = []

## 获取元素小数部分
dicMath = {}

for key in AtomicOccupationScoreDic:
    t = AtomicOccupationScoreDic[key]*AtomNumberof1
    dicMath[key] = modf(t)[0]

j = 0
k = 0
## 化整赋值
for key,value in dicMath.items():
    t = AtomicOccupationScoreDic[key]*AtomNumberof1

    if value >= 0.4:
        j = math.ceil(t)
        k = math.floor(AtomNumberofSingleAtom - j)
    else:
        j = math.floor(t)
        k = math.floor(AtomNumberofSingleAtom - j)
    AtomNumberof1Dic[' #'+key+f'-{sublattice1}'] = j
    AtomNumberof2Dic[' #'+key+f'-{sublattice2}'] = k


m1 = 0
m2 = 0
for value in AtomNumberof1Dic.values():
    m1 += value 
for value in AtomNumberof2Dic.values():
    m2 += value 

## 计算多余原子数
res1 = AtomNumberof1 - m1
res2 = AtomNumberof2 - m2

minKeyLst1=[]
maxKeyLst2=[]

# minValue=min(AtomNumberof1Dic.values())
minValue = min(filter(None, AtomNumberof1Dic.values())) ## 找出非零最小值
for i,j in AtomNumberof1Dic.items():
    if j==minValue:
        minKeyLst1.append(i)

maxValue=max(AtomNumberof2Dic.values())
for i,j in AtomNumberof2Dic.items():
    if j==maxValue:
        maxKeyLst2.append(i)


AtomNumberof1Dic[minKeyLst1[0]] = AtomNumberof1Dic[minKeyLst1[0]] + res1
AtomNumberof2Dic[maxKeyLst2[0]] = AtomNumberof2Dic[maxKeyLst2[0]]  + res2 


for key in AtomicOccupationScoreDic:
    for i in range(AtomNumberof1Dic[' #'+str(key)+f'-{sublattice1}']):
        AtomNumberof1List.append(' #'+str(key)+f'-{sublattice1}')
    for j in range(AtomNumberof2Dic[' #'+str(key)+f'-{sublattice2}']):
        AtomNumberof2list.append(' #'+str(key)+f'-{sublattice2}')

## 计算两种亚点阵的原子个数
def list_add(list1,list2):
    list3 = []
    for i in range(len(list1)):
        list3.append(list1[i]+list2[i])
    return list3

ls1 = []
ls2 = []

for i in AtomNumberof1Dic.values():
    ls1.append(i)

for i in AtomNumberof2Dic.values():
    ls2.append(i)

ls3 = []
ls3 = list_add(ls1,ls2)


## 1a
df1 = pd.read_excel(f'{workPath}\BCC-原子位置模板{strOcc}.xlsx',sheet_name=0)

df1 = df1.reindex(np.random.permutation(df1.index))
csv1 = df1.to_csv(f'{workPath}\{sublattice1}.csv',index=False)



df2 = pd.read_csv(f'{workPath}\\{sublattice1}.csv')


df2['元素'] = [i for i in AtomNumberof1List ]
df2.to_excel(f'{workPath}\\{sublattice1}.xlsx',index=False)

## 1b
df3 = pd.read_excel(f'{workPath}\BCC-原子位置模板{strOcc}.xlsx',sheet_name=1)

df3 = df3.reindex(np.random.permutation(df3.index))
csv2 = df3.to_csv(f'{workPath}\{sublattice2}.csv',index=False)

df4 = pd.read_csv(f'{workPath}\{sublattice2}.csv')
df4['元素'] = [i for i in AtomNumberof2list ]

df4.to_excel(f'{workPath}\\{sublattice2}.xlsx',index=False)

## 将1a和1b合并
df5 = pd.DataFrame()
df5 = pd.concat([df2,df4])


df5.sort_values(by='元素',inplace=True) ## 按距原子离排序

df5.to_excel(f'{workPath}\\{sublattice3}.xlsx',index=False) ## 输出不带索引的数据

## 获取元素名称字符串
lstName = [value for value in AtomicOccupationScoreDic.keys()]
strName = ''

for i in lstName:
    strName = strName + ''.join(i)
    

#####################################  POSCAR 文件生成  ##########################################

## poscar 配置信息
poscar1 = PosCar(f"{workPath}\\poscarTemplate\\{sublattice1}-{strOcc}-{AtomNumberof1}atoms-temp.vasp")
poscar2 = PosCar(f"{workPath}\\poscarTemplate\\{sublattice2}-{strOcc}-{AtomNumberof2}atoms-temp.vasp")
poscar3 = PosCar(f"{workPath}\\poscarTemplate\\{sublattice3}-{strOcc}-{atomSum}atoms-temp.vasp")


#  生产1a POSCAR文件

array1 = df2.to_numpy()

## poscar1 文件配置信息
poscar1.bases_const = 1.0
poscar1.bases = poscarBases
poscar1.atom_types = [value for value in AtomicOccupationScoreDic.keys()]
poscar1.atom_numbers = ls1
poscar1.data = array1


#  生产1b POSCAR文件

array2 = df4.to_numpy()

## poscar2 文件配置信息
poscar2.bases_const = 1.0
poscar2.bases = poscarBases
poscar2.atom_types = [value for value in AtomicOccupationScoreDic.keys()]
poscar2.atom_numbers = ls2
poscar2.data = array2


#  生产1a+1b POSCAR文件

array3 = df5.to_numpy()

## poscar3 文件配置信息
poscar3.bases_const = 1.0
poscar3.bases = poscarBases
poscar3.atom_types = [value for value in AtomicOccupationScoreDic.keys()]
poscar3.atom_numbers = ls3
poscar3.data = array3


def mkdir(path):
    isExists = os.path.exists(path)
    if not isExists:
        os.makedirs(path)
    else:
        pass
    
outPath = f"{workPath}\outPoscar"

mkdir(outPath)

## 输出文件
poscar1.tofile(f"{outPath}\\{strName}-{sublattice1}-{strOcc}-{AtomNumberof1}atoms.vasp")
poscar2.tofile(f"{outPath}\\{strName}-{sublattice2}-{strOcc}-{AtomNumberof2}atoms.vasp")
poscar3.tofile(f"{outPath}\\{strName}-{sublattice3}-{strOcc}-{atomSum}atoms.vasp")

print(f'{strName}-{strOcc}-计算完成！')

