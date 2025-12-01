from abc import abstractmethod
from ast import Index
from operator import index
from os import makedirs, name, path, remove
from typing import List, Tuple
import numpy as np
from numpy.lib import index_tricks
from numpy.lib.function_base import disp, percentile
import pandas as pd
from pandas.core.accessor import DirNamesMixin
from collections import Counter
from numpy.core.numeric import outer
import os

## step 0 配置文件

## 工作路径
excelFileName = '1'
pathWork = r'C:\\Users\\Karas\\Desktop\\1'

## 初始Excel表格路径
data = pd.read_excel(f'{pathWork}\\{excelFileName}.xlsx')

## 计算元素 Fe Mn Co Cr
atomList = ['Co','Cr','Fe']

## 原子距离
nearNumber = 2.99

## 输入占位数
# atomSeriesNumbers = [4,5,6,7,8,9,10]
atomSeriesNumbers = [6,7,8,9,10,11,12,13,14]

## 创建工作文件夹
def mkdir( path):
    # path = path.strip()
    # path = path.restrip("\\")

    isExists = os.path.exists(path)

    if not isExists:
        os.makedirs(path)
    else:
        pass


## datasets 输出文件路径
dataSets = f'{pathWork}\{excelFileName}-datasets'

## 输出原子距离文件路径
AtomDistance = f'{pathWork}\{excelFileName}-AtomDistance-1nn'

## 输出原子占位数文件路径
AtomDistanceSort = f'{pathWork}\{excelFileName}-AtomDistanceSort-1nn'

mkdir(dataSets)
makedirs(AtomDistance)
mkdir(AtomDistanceSort)


## Step 1  切分原子文件 ##

# 循环切割出每个原子的坐标文件
for ele in atomList:

    ## 读取数据
    fileOutName = f'{ele}Data.xlsx'
    pathOut = os.path.join(dataSets,fileOutName)

    ## 按原子类别切割
    df =  data.loc[data.iloc[:,-1].str.find(ele) == 1]
    df.to_excel(pathOut)

    print('step1 '  + ele +' over')
    

## Step 2  计算原子距离  ##

def atomDis(atom,dataSets,AtomDistance,nearNumber):
   
    ## 读取数据
    fileReadName = f'{atom}Data.xlsx'
    fileOutName = f'{atom}{atom}-Distance-1nn.xlsx'
    pathIn = os.path.join(dataSets,fileReadName)
    pathOut = os.path.join(AtomDistance,fileOutName)
    
    ## 读取数据
    df1 = pd.read_excel(pathIn)

    df3 = df1.loc[:,'x':'z']

    ## 转变成数组
    arrayA = df3.to_numpy()
    lenA = len(arrayA)

    outcomeArr = []
    k = 0  ## 计数器
    for i in range(lenA):
        k=k+1
        print(k)
        for j in range(i+1):
            dis = np.linalg.norm(arrayA[i]-arrayA[j])*112.5
            if i == j:
                pass
            else:
                if dis <= nearNumber: ## 原子距离判定条件
                    outcomeArr.append([arrayA[i][0],arrayA[i][1],arrayA[i][2],'#',(i+1),arrayA[j][0],arrayA[j][1],arrayA[j][2],'#',(j+1),dis])
                else:
                    pass
            
    df2 = pd.DataFrame(outcomeArr)
    df2.columns =[f'{atom}1_x',f'{atom}1_y',f'{atom}1_z','#',f'{atom}1',f'{atom}2_x',f'{atom}2_y',f'{atom}2_z','#',f'{atom}2','Distance']

    df2.sort_values(by='Distance',inplace=True) ## 按距原子离排序

    ## 输出文件
    df2.to_excel(pathOut)

# 循环计算同类原子之间的距离
for i in atomList:
    atomDis(i,dataSets,AtomDistance,nearNumber)
    print('step2 '+str(i)+' over')


## Step 3  计算原子配位数  ##


def atomDisSort(atom,atomSeriesNumbers,AtomDistance,AtomDistanceSort):
    
    ## 读取数据
    fileReadName = f'{atom}{atom}-Distance-1nn.xlsx'
    fileOutName = f'{atom}{atom}-SeriesNumbers-1nn-{atomSeriesNumbers}.csv'
    pathIn = os.path.join(AtomDistance,fileReadName)
    pathOut = os.path.join(AtomDistanceSort,fileOutName)
 
    df = pd.read_excel(pathIn)
    df2 = df.iloc[:,5:6] 
    # df2 = df.iloc[:,4:5]

    ## 转变成数组
    arrayA = df2.to_numpy()

    ## 转变为列表
    lst = []
    for i in arrayA:
        lst.append(int(i))

    ## 计算每个原子出现的个数，将原子序号和出现个数存入字典
    b = dict(Counter(lst))

    ## 计算原子配位数（连续数）符合要求的
    countdic = {key:value for key,value in b.items()if value == atomSeriesNumbers}

    ## 输出符合要求的key值（原子序号），并用列表存储
    outkeys = countdic.keys()
    outkeys_lst = []
    for i in outkeys:
        outkeys_lst.append(i)

    ## 判断列表是否为空
    if outkeys_lst == []:
        pass
    else:
        ## 用key值（原子序号）作为索引,回原始数据df中寻找符合要求的行，并添加到df5中
        df5 = pd.DataFrame()
        
        for j in range(len(outkeys_lst)):
            df4 = df.loc[df.loc[:,f'{atom}1'] == outkeys_lst[j]]
            df5 = df5.append(df4,ignore_index=True)

        # 输出文件
        df5.to_csv(pathOut)


# 遍历每个原子距离文件，找出符合指定配位数的原子对
for i in atomList:
    for j in atomSeriesNumbers:
        atomDisSort(i,j,AtomDistance,AtomDistanceSort)
        print('step 3 '+str(i)+'-'+str(j)+' over')


print('计算完成')