
##################################################################
##
## 该文件定义程序的介绍
## 小曲
## 2025/02/07
##
##################################################################



##################################################################

## 这个函数用来输出整个程序的介绍

def print_program_introduction():
    print(                                                                                                  )
    print( "------------------------------------------------------------------------------------------"     )  
    print(                                                                                                  )
    print( "     数值相对论计算程序 AMSS-NCKU  "                                                              )
    print( "     Numerical Relativity AMSS-NCKU  "                                                         )
    print(                                                                                                 )
    print( "     原程序作者：                     曹周键等 "                                                               )
    print( "     程序 Python 接口作者：小曲 "                                                                 )
    print(                                                                                                 )
    print( "     AMSS-NCKU 是一个数值相对论程序  "                                                            )
    print( "     本程序用来对双黑洞合并过程进行数值模拟，通过数值求解爱因斯坦方程得到双黑洞合并过程中引力场随时间的演化，" )
    print( " 从而得到黑洞的轨迹和释放引力波的信息。"                                                             )
    print(                                                                                                 )
    print( "     在数值方法上，本程序用使用有限差分方法对双黑洞合并过程进行数值模拟。程序中可以选择的差分方法有：2 阶 "   )
    print( " 差分、4 阶差分、6 阶差分、8 阶差分 "                                                               )
    print( "     程序中可以选择的微分方程为：BSSN 方程、Z4C 方程、BSSN 方程耦合标量场、BSSN 方程耦合电磁场 "          ) 
    print( "     程序中可以选择的网格类型为：方形网格、最外层带球壳的 shell patch 网格 "                            )
    print(                                                                                                  )
    print( "     除此之外，本程序还实现了 CPU 和 GPU 的混合运算"                                                 )
    print(                                                                                                 )
    print(                                                                                                 )
    print( "     Numerical Relativity AMSS-NCKU  "                                                                  )
    print(                                                                                                          )
    print( "     Author of AMSS-NCKU Code: Zhou-Jian Cao et al. "                                                   )
    print( "     Author of AMSS-NCKU Python Interface: Xiao Qu "                                                    )
    print(                                                                                                          )
    print( "     AMSS-NCKU is an open source numerical relativity code "                                            )
    print( "     It can be used to simulate the dynamical evolution on mergering process of black hole systems, "   )
    print( " calculating the variation of gravitational field, black holes' trajectories, and gravitational wave "  )
    print( " emissions through directly solving the Einstein field equations  "                                     )
    print(                                                                                                          )
    print( "     This AMSS-NCKU code uses the finite-difference method to evaluate the numerical simulation. The "  )
    print( " finite-difference schemes can be chosen as: 2nd order, 4th order, 6th order, 8th order. "              )
    print( "     The computation equation form in AMSS-NCKU code can be chosen as: BSSN equations, Z4C equations, " )
    print( " BSSN equations coupled with scalars (in f(R) theory), BSSN equations coupled with electromagnetic "    )
    print( " fields. "                                                                                              ) 
    print( "     The numerical grid system in this code includes: patch AMR grid, shell-patch AMR grid. "           )
    print(                                                                                                          )
    print( "     Furthermore, This code has fulfilled the CPU and GPU hybrid calculation. "                         )
    print(                                                                                                          )
    print( "------------------------------------------------------------------------------------------"             ) 
    print(                                                                                                          )
 
##################################################################

