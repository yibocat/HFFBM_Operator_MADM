import copy
import re
import networkx

import numpy as np
import pandas as pd

import mohusets.fuzzynumbers as fns
import mohusets.fuzzysets as fs
import mohusets.fuzzymeasure as fm

from matplotlib import pyplot as plt


def show_result(Alter, method, result, reverse=False):
    """
        展示决策结果
        输入：
            Alter 表示决策矩阵
            method （string 类型）表示方法名称
            result 表示由方法得到的计算结果
            reverse 表示排序方式，False 表示升序，True 表示降序
    """
    suppliers = []
    for i in range(len(Alter)):
        suppliers.append('A' + str(i + 1))
    dic = dict(zip(suppliers, result))
    rank_table = pd.DataFrame(columns=suppliers)
    rank_table.loc[method] = (list(dic.values()))
    sort = []
    for a in range(len(suppliers)):
        sort.append(sorted(dic.items(), key=lambda d: d[1], reverse=reverse)[a][0])
    m = '%s' % sort[0]
    for aa in sort[1:]:
        m = m + '>' + aa
    rank_table['Ranking Order'] = m
    return rank_table


# 最大偏差法求权重
def maxdev(alter, lam=1, tao=1, indeterminacy=True):
    A = alter.T
    distance_list = [0 for i in range(len(alter[0]))]  # 计算的距离矩阵
    for att in range(len(alter[0])):
        sum_dis = []
        for i, d1 in enumerate(A[att]):
            for j, d2 in enumerate(A[att][:i]):
                sum_dis.append(fns.generalized_distance(d1, d2, lam, tao, indeterminacy))
        sum_dis = np.asarray(sum_dis)
        distance_list[att] = sum_dis.sum()  # 计算每个属性下的偏差和
    distance_list = np.asarray(distance_list)  # 转换为 ndarray 类型
    sum_distance_list = distance_list.sum()

    weight_list = []
    for att in distance_list:
        weight_list.append(att / sum_distance_list)
    return np.asarray(weight_list)


def ranking_dm(operator, A, wei, reverse=False, **params):
    """
    将费马模糊决策矩阵通过有序加权聚合算子进行有序转换，转换成有序费马模糊决策矩阵
    :param operator:    函数形参，传递函数，表示使用哪种聚合算子，分别为 ffhhwa 和 ffhhwg
    :param reverse:     排序方式：True 表示降序，False 表示升序, 默认为降序
    :param A:           费马决策矩阵
    :param wei:         属性权重
        **param 表示形参，主要是 opera 函数的参数，用 “a=1,b=2,c=3” 表示，具体参考聚合算子的参数
    :return:            返回使用聚合算子排序后的费马模糊决策矩阵
    """

    def _ranking(opera, alter, w, rev, **param):
        """
        作为有序加权，对权重和 ffn 进行扩展，再按照聚合算子的得分值进行排序
        相当于有序加权算子
        :param rev:     排序方式：True 表示降序，False 表示升序
        :param alter:   待聚合的费马模糊数集合
        :param w:       权重，这里是属性权重
        :return: 排序后的费马模糊数集合
        """
        # 扩展权重
        sw = []
        for i in range(len(w)):
            ww = [w[i] for j in range(len(w))]
            sw.append(ww)

        # 扩展聚合 ffn 列表
        al = []
        for i in range(len(alter)):
            ll = [alter[i] for j in range(len(alter))]
            al.append(ll)

        # 计算聚合得分
        score_list = []
        for att in range(len(al)):
            # print(ffhhwa(al[att],sw[att],lam))            # 打印聚合后的费马模糊数
            score_list.append(opera(al[att], sw[att], **param).score)

        # 排序
        sorted_list = sorted(dict(zip(alter, score_list)).items(), key=lambda d: d[1], reverse=rev)
        return sorted_list

    sort_dm = []
    for i in range(len(A)):
        sort_dm.append(_ranking(operator, A[i], wei, reverse, **params))

    so = []
    for i in range(len(sort_dm)):
        st = []
        for j in range(len(sort_dm[i])):
            st.append(sort_dm[i][j][0])
        so.append(st)
    # return np.asarray(so)
    return fs.asfuzzyset(so)

