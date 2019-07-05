# Allocation
goods allocation

transportation.py     some functions


allocation2.py          handle the allocation problem
                        using command: python allocation2.py surplus_PATH lack_PATH costs_PATH results_PATH solve_method iter_num
                        argv:   1, surplus_PATH: 供给量csv文件路径，第一行和第一列分别是行列index，剩余是N×K的矩阵，第i行第j列元素是i号店铺可以提供的j
                                号sku的数量； 
                                2, lack_PATH： 需求量csv文件路径，格式同上 
                                3, costs_PATH： 运输时间成本csv文件路径，第一行和第一列分别是行列index，剩余是N×N的矩阵，第i行j列元素是i号店铺到j号店铺的运输时间；
                                4, results_PATH： 输出文件路径，输出csv文件，第一行是列index，第一列是发货店铺编号(i)，第二列是收货店铺编号(j)，第k+2列是第k种sku的运输量，如第m行n列(n>=2)表示从 result[m][0] 号店铺向 result[m][1] 号店铺运输 result[m][n]件 n-2 种sku
                                注意的是该文件没有保存行index
                                5, solve_method： 求解方法，建议选1 
                                6, iter_num： 若选择方法1，则可以设置最大迭代次数，整数

                                其中1,2,3,5,6是输入，4是输出


上述函数需要的输入文件可以通过data_ps.py生成

data_ps.py              第7行： supply = pd.read_csv('./data/supply_0623.csv', encoding='gbk') 需要修改输入文件，是公司提供的原始文件
                        第9行： pre_results = pd.read_csv('./data/623_results_nm.csv', encoding='utf-8')  需要改成李刚提供的结果
                        第69行：store_Frame.to_csv('~/production/data_623/store_list_nm.csv')  需要修改输出文件名，后续72,75,78,81行同样