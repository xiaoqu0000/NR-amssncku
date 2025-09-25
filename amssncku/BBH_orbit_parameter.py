
##############################################################################################

## 这个程序用于设定双黑洞旋转时的参数
## 小曲
## 2024/04/05
## 修改 2025/02/10

## 可以用于 AMSS-NCKU 或 Einstein Toolkit 程序的输入

##############################################################################################


import AMSS_NCKU_Input as input_data
import math
import os
import sympy
import numpy      
import derivative    ## 数值求导


##############################################################################################

## 根据输入文件设置每个黑洞的真实角动量 

angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )   ## 初始化每个黑洞的自旋角动量

for i in range(input_data.puncture_number):
    if ( input_data.Symmetry == "equatorial-symmetry" ):
        angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
    elif ( input_data.Symmetry == "no-symmetry" ):
        angular_momentum_BH[i] = (input_data.parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i] 

## 设置两黑洞质量
## 为了跟文献的记号一致，这里要求  M1 >= M2

M1 = input_data.parameter_BH[0,0]
M2 = input_data.parameter_BH[1,0] 

## 设置两黑洞的无量纲自旋

S1 = angular_momentum_BH[0] / M1**2
S2 = angular_momentum_BH[1] / M2**2

## 设置两黑洞质心系中轨道半长轴和偏心率

D0 = input_data.Distance
e0 = input_data.e0

##############################################################################################


##############################################################################################

## 该函数设定了双黑洞圆轨道旋转时的轨道参数
## 作者：小曲
## 使用后牛顿近似获得准圆轨道
## 更新到 3 阶后牛顿近似

## 自变量
## 双黑洞质量 M1 和 M2 （注意一定要 M1 > M2）
## 双黑洞自旋 S1 和 S2 （它们必须为 numpy 的 3 个元素构成的向量，否则无法用 numpy.dot ）
## 初始间距 D0
## 轨道椭率 e0 （现在只包含了圆轨道的后牛顿公式，以后可以加入带椭率的公式）

def generate_BBH_orbit_parameters( M1, M2, S1, S2, D0, e0 ):

    print(                                                )
    print(                                                )
    print( " 根据双星系统有效单体模型计算出双星轨道的特征量 " )
    print(                                                ) 

    print(                                                                    )
    print( f" 已输入双星的质量为：       M1 = {M1}  M2 = {M2} "                 )
    print( f" 已输入双星的无量纲自旋为： S1 = {S1}  S2 = {S2} "                 )
    print( f" 已输入双星轨道半长轴和偏心率为： a0 = D0/2 = {D0/2.0}  e0 = {e0} " )
    print(                                                                    )
    print(  " 下面开始计算 "                                                   )
    print(                                                                    )

    ##################################################

    ## 求出双星轨道的质量比，约化质量，无量纲质量等
    M_total = M1 + M2
    print( f" 双星质量： M1 = {M1}  M2 = {M2}  牛顿力学总质量：M_total = {M_total} " )

    ## 求出约化质量
    M_mu = M1 * M2 / M_total
    print( " 单体轨道约化质量为：M_mu = ", M_mu )

    ## 设置质量比
    ## 特别注意，为了跟 TwoPuncture 的计算一致，这里要求  M1 >= M2
    m_q   = M1 / M2
    m_eta = m_q / ( (1.0 + m_q)**2 )
    print( " 无量纲质量比为：Q = M1 / M2 = ", m_q )

    ## 设置无量纲质量
    m1   = M1 / M_total
    m2   = M2 / M_total
    m_mu = M_mu / M_total
    print( f" 无量纲质量 ：m1 ={m1}  m2 = {m2}  m_mu = {m_mu} " )
    print(  " 无量纲约化质量为： m_eta = Q / (1+Q)^2 = ", m_eta  )

    ##################################################    
    
    ## 由质心系半长轴和偏心率求出 t = 0 时双星对应的轨道参数

    ## 根据可以经典力学，双星的轨道可以等效为约化质量在中心力场的运动
    ## 单体轨道半径 r，相位角 phi，半长轴 a，偏心率 e 存在关系
    ## r = a*(1-e^2)/(1+e*cos(phi))
    ## 相位角 phi = 0 或 phi = pi 对应半长轴和半短轴

    ## 求出轨道半长轴与半短轴
    a0 = D0 / 2.0
    a_long  = a0 * (1.0 + e0)
    a_short = a0 * (1.0 - e0)

    ## 根据质心系本身的定义可以求出双星各自的半长轴
    ## M1/M2 = a2/a1
    R10 = D0 * M2 / M_total
    R20 = D0 * M1 / M_total

    print(                        )
    print( " 求出双星轨道坐标参数 " )
    print(                        )
    print( " 因为坐标轴可以任意选，不妨设 t = 0 时刻 y 轴方向正好与轨道半长轴重合，双星共同沿着 z 方向转动 " )
    print( " 此时坐标 y 关联径向坐标 R，坐标 x 关联角向坐标 phi，z 坐标在 t = 0 为零 "                    )
    print( " 同样，t = 0 时，动量 Py 关联径向坐标 Pr，动量 Px 关联角向动量 P_phi，动量 Pz 为零 "           )
    ## print( " 请注意，这并不意味着 z 坐标在演化时间很久后永远仍然为零，极端情况下（比如自旋角动量方向不在 z 轴方向）轨道可能偏离 xy 平面 " ) 
    print()

    ## 初始化位置坐标
    position1 = [0.0, 0.0, 0.0]
    position2 = [0.0, 0.0, 0.0]

    ## 初始化动量坐标
    momentum1 = [0.0, 0.0, 0.0]
    momentum2 = [0.0, 0.0, 0.0]


    position1[1] =   R10    
    position2[1] = - R20

    print( f" 双星在 t=0 时刻的坐标为："            )
    print( f" Y1 = {position1}  Y2 = {position2}" )
    print(                                        )


    ########################################

    ## 下面求出双星轨道的旋转角速度

    R = D0
    epsilon = (m1 + m2) / R 
    print( " 后牛顿展开参数 epsilon =  M/r = ", epsilon )
    print(                                             )
    ## 双星系统引力波具有标度不变性，可用无量纲总质量 m = m1 + m2 = 1 进行计算，且这里的初始间距 D0 默认也是以总质量为单位（实际上给定的是 真实距离/M_total）

    ## 3 PN 后牛顿近似的结果
    ## 注意：以下公式对应圆轨道，有偏心率的情况还没有考虑进来
    ## 根据文献  
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower 
    ## Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity, 
    ## arXiv:1702.00872[gr-qc]
    ## 还要特别注意，文献 arXiv:1702.00872[gr-qc] 中的质量比 q 是小于 1 的，因此文献中的 S1 对应这里的 S2 

    ## 设定双星旋转角速度（以后可以修改这个地方）

    Omega_0PN = epsilon**1.5

    Omega_correction_1PN  = - 0.5 * epsilon * ( ( 3.0*(m_q**2) + 5.0*m_q + 3.0 ) / (1.0+m_q)**2 ) 
    Omega_correction_15PN = - 0.25 * ( epsilon**(1.5) )                               \
                              * (   (3.0 + 4.0*m_q) * m_q * S2[2] / ( (1.0+m_q)**2 )  \
                                  + (3.0*m_q + 4.0)       * S1[2] / ( (1.0+m_q)**2 )        \
                                )
    Omega_correction_2PN  = epsilon**2 \
                            * (   (1.0/16.0) * ( 24.0*(m_q**4) + 103.0*(m_q**3) + 164.0*(m_q**2) + 103.0*m_q + 24.0 ) \
                                   / ( (1.0+m_q)**4 )                                 \
                                - 1.5  * (m_q**2) * (S2[0]**2) / ( (1.0+m_q)**2 )     \
                                + 0.75 * (m_q**2) * (S2[1]**2) / ( (1.0+m_q)**2 )     \
                                + 0.75 * (m_q**2) * (S2[2]**2) / ( (1.0+m_q)**2 )     \
                                - 3.0  * m_q * S1[0] * S2[0] / ( (1.0+m_q)**2 )       \
                                + 1.5  * m_q * S1[1] * S2[1] / ( (1.0+m_q)**2 )       \
                                + 1.5  * m_q * S1[2] * S2[2] / ( (1.0+m_q)**2 )       \
                                - 1.5  * (S1[0]**2) / ( (1.0+m_q)**2 )                \
                                + 0.75 * (S1[1]**2) / ( (1.0+m_q)**2 )                \
                                + 0.75 * (S1[2]**2) / ( (1.0+m_q)**2 )
                              )
    Omega_correction_25PN = (3.0/16.0) * (epsilon**(2.5))   \
                            * (   S2[2] * m_q * ( 16.0*(m_q**3) + 30.0*(m_q**2) + 34.0*m_q + 13.0 ) / ( (1.0+m_q)**4 )  \
                                + S1[2]       * ( 13.0*(m_q**3) + 34.0*(m_q**2) + 30.0*m_q + 16.0 ) / ( (1.0+m_q)**2 )  \
                              )
    Omega_correction_3PN  = epsilon**3 \
                            * (   (167.0/128.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )           \
                                - (   120.0*(m_q**6) + 2744.0*(m_q**5) + 10049.0*(m_q**4)         \
                                    + 14820.0*(m_q**3) + 10049.0*(m_q**2) + 2744.0*m_q + 120.0    \
                                   ) / ( 96.0 * ((1.0+m_q)**6) )                                  \
                                + (1.0/16.0) * (m_q**2) * (S2[0]**2) * ( 76.0*(m_q**2) + 180.0*m_q + 155.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/8.0)  * (m_q**2) * (S2[1]**2) * ( 43.0*(m_q**2) + 85.0*m_q  + 55.0 )    / ( (1.0+m_q)**4 ) \
                                - (1.0/32.0) * (m_q**2) * (S2[2]**2) * ( 2.0*m_q + 5.0 ) * ( 14.0*m_q + 27.0 ) / ( (1.0+m_q)**4 ) \
                                + (1.0/16.0)            * (S1[0]**2) * ( 155.0*(m_q**2) + 180.0*m_q + 76.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/8.0)             * (S1[1]**2) * ( 55.0*(m_q**2)  + 85.0*m_q  + 43.0 )   / ( (1.0+m_q)**4 ) \
                                - (1.0/32.0)            * (S1[2]**2) * ( 27.0*m_q + 14.0 ) * ( 5.0*m_q + 2.0 ) / ( (1.0+m_q)**4 ) \
                                + (1.0/8.0)  * m_q * S1[0] * S2[0] * ( 120.0*(m_q**2) + 187.0*m_q + 120.0 ) / ( (1.0+m_q)**4 )    \
                                - 0.25       * m_q * S1[1] * S2[1] * ( 54.0*(m_q**2) + 95.0*m_q   + 54.0 )  / ( (1.0+m_q)**4 )    \
                                - (1.0/16.0) * m_q * S1[2] * S2[2] * ( 96.0*(m_q**2) + 127.0*m_q  + 96.0 )  / ( (1.0+m_q)**4 )    \
                              )

    Omega_1PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN )
    Omega_15PN = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN )
    Omega_2PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN                            \
                             )
    Omega_25PN = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN + Omega_correction_25PN    \
                             )
    Omega_3PN  = Omega_0PN * ( 1.0 + Omega_correction_1PN + Omega_correction_15PN    \
                                   + Omega_correction_2PN + Omega_correction_25PN    \
                                   + Omega_correction_3PN                            \
                             )

    print(                                                   )
    print( " 0   阶后牛顿近似下双星旋转的角速度 = ", Omega_0PN  )
    print( " 1   阶后牛顿近似下双星旋转的角速度 = ", Omega_1PN  )
    print( " 1.5 阶后牛顿近似下双星旋转的角速度 = ", Omega_15PN )
    print( " 2   阶后牛顿近似下双星旋转的角速度 = ", Omega_2PN  )
    print( " 2.5 阶后牛顿近似下双星旋转的角速度 = ", Omega_25PN )
    print( " 3   阶后牛顿近似下双星旋转的角速度 = ", Omega_3PN  )
    print(                                                   )

    ########################################

    ## 设定双星旋转的角向动量（以后可以修改这个地方）

    Pt_0PN = (epsilon**0.5) * m_q / ( (1+m_q)**2.0 ) 

    Pt_correction_1PN  = 2.0 * epsilon
    Pt_correction_15PN = epsilon**1.5 \
                         * ( - 0.75 * (3.0 + 4.0*m_q) * m_q * S2[2] / ( (1.0+m_q)**2 ) \
                             - 0.75 * (3.0*m_q + 4.0)       * S1[2] / ( (1.0+m_q)**2 )       \
                           )
    Pt_correction_2PN  = epsilon**2.0 \
                         * (  (1.0/16.0) * ( 42.0*(m_q**2) + 41.0*m_q + 42.0 ) / ( (1.0+m_q)**2 )   \
                             - 1.5  * (m_q**2) * (S2[0]**2) / ( (1.0+m_q)**2 )                      \
                             + 0.75 * (m_q**2) * (S2[1]**2) / ( (1.0+m_q)**2 )                      \
                             + 0.75 * (m_q**2) * (S2[2]**2) / ( (1.0+m_q)**2 )                      \
                             - 3.0  * m_q * S1[0] * S2[0] / ( (1.0+m_q)**2 )                        \
                             + 1.5  * m_q * S1[1] * S2[1] / ( (1.0+m_q)**2 )                        \
                             + 1.5  * m_q * S1[2] * S2[2] / ( (1.0+m_q)**2 )                        \
                             - 1.5  * (S1[0]**2) / ( (1.0+m_q)**2 )                                 \
                             + 0.75 * (S1[1]**2) / ( (1.0+m_q)**2 )                                 \
                             + 0.75 * (S1[2]**2) / ( (1.0+m_q)**2 )                                 \
                           )
    Pt_correction_25PN = epsilon**2.5 \
                         * ( - (1.0/16.0) * ( 72.0*(m_q**3) + 116.0*(m_q**2) + 60.0*m_q + 13.0 ) \
                                          * m_q * S2[2] / ( (1.0+m_q)**4 )                       \
                             - (1.0/16.0) * ( 13.0*(m_q**3) + 60.0*(m_q**2) + 116.0*m_q + 72.0 ) \
                                          * S1[2] / ( (1.0+m_q)**4 )                             \
                           )
    Pt_correction_3PN  = epsilon**3.0 \
                         * (   (163.0/128.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )                                   \
                             + (1.0/32.0) * ( 120.0*(m_q**4) - 659.0*(m_q**3) - 1532.0*(m_q**2) - 659.0*m_q + 120.0 )  \
                                / ( (1.0+m_q)**4 )                                                                     \
                             - (1.0/16.0) * (S2[0]**2) * (m_q**2) * ( 80.0*(m_q**2)             - 59.0 ) / ( (1.0+m_q)**4 ) \
                             - 0.5        * (S2[1]**2) * (m_q**2) * (        m_q**2  + 10.0*m_q + 8.0  ) / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * (S2[2]**2) * (m_q**2) * ( 128.0*(m_q**2) + 56.0*m_q - 27.0 ) / ( (1.0+m_q)**4 ) \
                             - (1.0/16.0) * (S1[0]**2)            * ( 80.0 - 59.0*(m_q**2) )             / ( (1.0+m_q)**4 ) \
                             - 0.5        * (S1[1]**2)            * ( 8.0*(m_q**2) + 10.0*m_q + 1.0 )    / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * (S1[2]**2)            * ( 128.0 + 56.0*m_q - 27.0*(m_q**2) ) / ( (1.0+m_q)**4 ) \
                             + (1.0/8.0)  * S1[0] * S2[0] * m_q   * ( 12.0*(m_q**2) + 35.0*m_q + 12.0 )  / ( (1.0+m_q)**4 ) \
                             - 0.25       * S1[1] * S2[1] * m_q   * ( 27.0*(m_q**2) + 58.0*m_q + 27.0 )  / ( (1.0+m_q)**4 ) \
                             + (1.0/32.0) * S1[2] * S2[2] * m_q   * ( 60.0*(m_q**2) + 13.0*m_q + 60.0 )  / ( (1.0+m_q)**4 ) \
                           )

    Pt_1PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN )
    Pt_15PN = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN )
    Pt_2PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN                         \
                       )
    Pt_25PN = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN + Pt_correction_25PN    \
                       )
    Pt_3PN  = Pt_0PN * ( 1.0 + Pt_correction_1PN + Pt_correction_15PN    \
                             + Pt_correction_2PN + Pt_correction_25PN    \
                             + Pt_correction_3PN
                       )

    print(                                                  )
    print( " 0   阶后牛顿近似下双星旋转的角向动量 = ", Pt_0PN  )
    print( " 1   阶后牛顿近似下双星旋转的角向动量 = ", Pt_1PN  )
    print( " 1.5 阶后牛顿近似下双星旋转的角向动量 = ", Pt_15PN )
    print( " 2   阶后牛顿近似下双星旋转的角向动量 = ", Pt_2PN  )
    print( " 2.5 阶后牛顿近似下双星旋转的角向动量 = ", Pt_25PN )
    print( " 3   阶后牛顿近似下双星旋转的角向动量 = ", Pt_3PN  )
    print(                                                  )

    ########################################

    ## 下面求出双星系统的 ADM 质量
    ## 根据文献  
    ## Antoni Ramos-Buades, Sascha Husa, and Geraint Pratten
    ## Simple procedures to reduce eccentricity of binary black hole simulations 
    ## arXiv:1810.00036[gr-qc]

    ############################

    ## 定义 ADM 质量函数
    ## 该函数用 M/R 来展开

    def M_ADM(r):
        
        mass     =  m1 + m2
        epsilon0 = (m1 + m2) / r
        
        adm_correction_0PN  = - 0.5 * epsilon0 * m_q / ( (1+m_q)**2 ) 
        adm_correction_1PN  =  (1.0/8.0) * ( epsilon0**2 )                                   \
                                * m_q * ( 7.0*(m_q**2) + 13.0*m_q + 7.0 ) / ( (1.0+m_q)**4 ) 
        adm_correction_15PN = - 0.25 * ( epsilon0**2.5 )                                     \
                                * (   m_q**2 * ( 3.0 + 4.0*m_q ) * S2[2] / ( (1.0+m_q)**4 )  \
                                    + m_q    * ( 3.0*m_q + 4.0 ) * S1[2] / ( (1.0+m_q)**4 )  \
                                  )
        adm_correction_2PN  = epsilon0**3  \
                              * (   (1.0/16.0) * m_q * ( 9.0*m_q**4 + 16.0*(m_q**3) + 13.0*(m_q**2) + 16.0*m_q + 9.0 ) \
                                     / ( (1.0+m_q)**6 )                                                                \
                                  - 0.5  * (S2[0]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  + 0.25 * (S2[1]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  + 0.25 * (S2[2]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                  - 1.0  * S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  + 0.5  * S1[1] * S2[1] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  + 0.5  * S1[2] * S2[2] * (m_q**2) / ( (1.0+m_q)**4 )   \
                                  - 0.5  * (S1[0]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  + 0.25 * (S1[1]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  + 0.25 * (S1[2]**2) * m_q / ( (1.0+m_q)**4 )           \
                                )
        adm_correction_25PN = - (1.0/16.0) * epsilon0**3.5  \
                                * (   S2[2] * (m_q**2) * ( 32.0*(m_q**3) + 42.0*(m_q**2) + 14.0*m_q +1.0 ) / ( (1.0+m_q)**6 )  \
                                    + S1[2] * m_q      * ( m_q**3 + 14.0*(m_q**2) + 42.0*m_q + 32.0 )      / ( (1.0+m_q)**6 )  \
                                  )
        adm_correction_3PN  = epsilon0**4  \
                              * (   (81.0/128.0) * (math.pi**2) * (m_q**2) / ( (1.0+m_q)**4 )    \
                                  + (     537.0*(m_q**6) - 3497.0*(m_q**5) - 18707.0*(m_q**4)    \
                                      - 29361.0*(m_q**3) - 18707.0*(m_q**2) - 3497.0*m_q + 537.0 \
                                    ) * (m_q/384.0) / ( (1.0+m_q)**8 )                           \
                                  - (1.0/16.0) * (S2[0]**2)    * (m_q**3) * ( 52.0*(m_q**2) + 12.0*m_q - 25.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/8.0)  * (S2[1]**2)    * (m_q**3) * (       m_q**2  - 17.0*m_q - 15.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S2[2]**2)    * (m_q**3) * ( 50.0*(m_q**2) + 38.0*m_q + 3.0  ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S1[0]**2)    * m_q      * ( 25.0*(m_q**2) - 12.0*m_q - 52.0 ) / ( (1.0+m_q)**6 )  \
                                  - (1.0/8.0)  * (S1[1]**2)    * m_q      * ( 15.0*(m_q**2) + 17.0*m_q -  1.0 ) / ( (1.0+m_q)**6 )  \
                                  + (1.0/16.0) * (S1[2]**2)    * m_q      * (  3.0*(m_q**2) + 38.0*m_q + 50.0 ) / ( (1.0+m_q)**6 )  \
                                  + (9.0/8.0)  * S1[0] * S2[0] * (m_q**3)                                       / ( (1.0+m_q)**6 )  \
                                  - (3.0/4.0)  * S1[1] * S2[1] * (m_q**2) * (  4.0*(m_q**2) +  9.0*m_q + 4.0)   / ( (1.0+m_q)**6 )  \
                                  + (3.0/8.0)  * S1[2] * S2[2] * (m_q**2) * ( 10.0*(m_q**2) + 21.0*m_q + 10.0 ) / ( (1.0+m_q)**6 )  \
                                )

        ADM_0PN  = mass * ( 1.0 + adm_correction_0PN )
        ADM_1PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN )
        ADM_15PN = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN 
                          )
        ADM_2PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                          )
        ADM_25PN = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                                + adm_correction_25PN                           \
                          )
        ADM_3PN  = mass * ( 1.0 + adm_correction_0PN  + adm_correction_1PN      \
                                + adm_correction_15PN + adm_correction_2PN      \
                                + adm_correction_25PN + adm_correction_3PN      \
                          )

        return ADM_0PN, ADM_1PN, ADM_15PN, ADM_2PN, ADM_25PN, ADM_3PN
    
    ############################

    ## 定义另一个 ADM 质量函数
    ## 该函数用旋转角速度 Omega 来展开
    ## 根据文献
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower 
    ## Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity, 
    ## arXiv:1702.00872[gr-qc]
    ## 还要特别注意，文献 arXiv:1702.00872[gr-qc] 中的质量比 q 是小于 1 的，因此文献中的 S1 对应这里的 S2 

    def M_ADM_another(Omega):
        
        mass     =  m1 + m2
        epsilon0 = mass * Omega
        
        adm_correction_0PN  = 0.0
        adm_correction_1PN  = - 0.5 * ( epsilon0**(2.0/3.0) )
        adm_correction_15PN = ( epsilon0**(4.0/3.0) )   \
                              * (1.0/24.0) * ( 9.0*(m_q**2) + 19.0*m_q + 9.0 ) / ( (1.0+m_q)**2 )
        adm_correction_2PN  = - (1.0/3.0) * ( epsilon0**(5.0/3.0) )   \
                                          * (   S2[2] * m_q * ( 4.0*m_q + 3.0 ) / ( (1.0+m_q)**2 ) \
                                              + S1[2] *       ( 3.0*m_q + 4.0 ) / ( (1.0+m_q)**2 ) \
                                            )                                                      \
                              + ( epsilon0**2 )   \
                                * (   (1.0/48.0) * (   81.0*(m_q**4) + 267.0*(m_q**3)     \
                                                     + 373.0*(m_q**2) + 267.0*m_q + 81.0  \
                                                   ) / ( (1.0+m_q)**4 )                   \
                                    -       (S2[0]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    + 0.5 * (S2[1]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    + 0.5 * (S2[2]**2) * (m_q**2) / ( (1.0+m_q)**2 )      \
                                    -       (S1[0]**2) / ( (1.0+m_q)**2 )                 \
                                    + 0.5 * (S1[1]**2) / ( (1.0+m_q)**2 )                 \
                                    + 0.5 * (S1[2]**2) / ( (1.0+m_q)**2 )                 \
                                    - 2.0 * S1[0] * S2[0] * m_q / ( (1.0+m_q)**2 )        \
                                    +       S1[1] * S2[1] * m_q / ( (1.0+m_q)**2 )        \
                                    +       S1[2] * S2[2] * m_q / ( (1.0+m_q)**2 )        \
                                  )
        adm_correction_25PN = - (1.0/18.0) * ( epsilon0**(7.0/3.0) )   \
                                * (   S2[2] * m_q * ( 72.0*(m_q**3) + 140.0*(m_q**2) + 96.0*m_q + 27.0 ) / ( (1.0+m_q)**4 ) \
                                    + S1[2]       * ( 27.0*(m_q**3) + 96.0*(m_q**2) + 140.0*m_q + 72.0 ) / ( (1.0+m_q)**4 ) \
                                  )
        adm_correction_3PN  = ( epsilon0**(8.0/3.0) )   \
                              * (   (205.0/192.0) * (math.pi**2) * m_q / ( (1.0+m_q)**2 )              \
                                  + (    54675.0*(m_q**6) + 18045.0*(m_q**5) - 411525.0*(m_q**4)       \
                                      - 749755.0*(m_q**3) - 411525.0*(m_q**2) + 18045.0*m_q + 54675.0  \
                                    ) / ( 10368.0 * ( (1.0+m_q)**6 ) )                                 \
                                  - (5.0/24.0) * (S2[0]**2) * (m_q**2) * ( 20.0*(m_q**2) + 4.0*m_q - 11.0 ) / ( (1.0+m_q)**4 )  \
                                  - (5.0/12.0) * (S2[1]**2) * (m_q**2) * (       m_q**2  + 9.0*m_q + 7.0  ) / ( (1.0+m_q)**4 )  \
                                  + (5.0/36.0) * (S2[2]**2) * (m_q**2) * ( 13.0*(m_q**2) - 3.0*m_q - 9.0  ) / ( (1.0+m_q)**4 )  \
                                  + (5.0/4.0)  * S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )                                    \
                                  - (6.0/5.0)  * S1[1] * S2[1] * m_q * ( 2.0*m_q + 3.0 ) * ( 3.0*m_q + 2.0 ) / ( (1.0+m_q)**4 ) \
                                  + (5.0/18.0) * S1[2] * S2[2] * m_q * ( 3.0*(m_q**2) + 7.0*m_q + 3.0 ) / ( (1.0+m_q)**4 )      \
                                  + (1.0/24.0) * (S1[0]**2) * ( 55.0*(m_q**2) - 20.0*m_q - 100.0 ) / ( (1.0+m_q)**4 )           \
                                  - (1.0/12.0) * (S1[1]**2) * ( 35.0*(m_q**2) + 45.0*m_q + 5.0   ) / ( (1.0+m_q)**4 )           \
                                  - (1.0/36.0) * (S1[2]**2) * ( 45.0*(m_q**2) + 15.0*m_q - 65.0  ) / ( (1.0+m_q)**4 )           \
                                )
        
        ADM_0PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * adm_correction_0PN )
        ADM_1PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                              ) 
                          )
        ADM_15PN = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN   \
                                                               + adm_correction_1PN   \
                                                               + adm_correction_15PN  \
                                                              ) 
                          )
        ADM_2PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                              )
                          )
        ADM_25PN = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                               + adm_correction_25PN \
                                                              )
                          )
        ADM_3PN  = mass * ( 1.0 + ( m_q / ( (1+m_q)**2 ) ) * (   adm_correction_0PN  \
                                                               + adm_correction_1PN  \
                                                               + adm_correction_15PN \
                                                               + adm_correction_2PN  \
                                                               + adm_correction_25PN \
                                                               + adm_correction_3PN  \
                                                              )
                          )

        return ADM_0PN, ADM_1PN, ADM_15PN, ADM_2PN, ADM_25PN, ADM_3PN
    
    ############################
    
    ## 定义 ADM 质量函数对 r 的导数

    def dADM_dr(r):
        
        mass     =  m1 + m2
        epsilon0 = (m1 + m2) / r
        
        dADM_correction_0PN  = ( - (epsilon0**2) / mass ) * ( - 0.5 * m_q / ((1+m_q)**2) ) 
        dADM_correction_1PN  = ( - 2.0 * (epsilon0**3) / mass )                                           \
                                * (1.0/8.0) * m_q * ( 7.0*(m_q**2) + 13.0*m_q + 7.0 ) / ( (1.0+m_q)**4 ) 
        dADM_correction_15PN = ( - 2.5 * (epsilon0**3.5) / mass )                            \
                                * ( - 0.25 )                                                 \
                                * (   m_q**2 * ( 3.0 + 4.0*m_q ) * S2[2] / ( (1.0+m_q)**4 )  \
                                    + m_q    * ( 3.0*m_q + 4.0 ) * S1[2] / ( (1.0+m_q)**4 )  \
                                  )
        dADM_correction_2PN  = ( - 3.0 * (epsilon0**4) / mass )   \
                                * (   (1.0/16.0) * m_q * ( 9.0*(m_q**4) + 16.0*(m_q**3) + 13.0*(m_q**2) + 16.0*m_q + 9.0 ) \
                                      / ( (1.0+m_q)**6 )                                                                   \
                                    - 0.5  * (S2[0]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    + 0.25 * (S2[1]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    + 0.25 * (S2[2]**2) * (m_q**3) / ( (1.0+m_q)**4 )      \
                                    - 1.0  *  S1[0] * S2[0] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    + 0.5  *  S1[1] * S2[1] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    + 0.5  *  S1[2] * S2[2] * (m_q**2) / ( (1.0+m_q)**4 )  \
                                    - 0.5  * (S1[0]**2) * m_q / ( (1.0+m_q)**4 )           \
                                    + 0.25 * (S1[1]**2) * m_q / ( (1.0+m_q)**4 )           \
                                    + 0.25 * (S1[2]**2) * m_q / ( (1.0+m_q)**4 )           \
                                  )
        dADM_correction_25PN = ( - 3.5 * (epsilon0**4.5) / mass )  \
                                * ( - 1.0/16.0 )                   \
                                * (   S2[2] * (m_q**2) * ( 32.0*(m_q**3) + 42.0*(m_q**2) + 14.0*m_q + 1.0  ) / ( (1.0+m_q)**6 )  \
                                    + S1[2] * m_q      * (       m_q**3  + 14.0*(m_q**2) + 42.0*m_q + 32.0 ) / ( (1.0+m_q)**6 )  \
                                  )
        dADM_correction_3PN  = ( - 4.0 * (epsilon0**5) / mass )   \
                                * (   (81.0/128.0) * (math.pi**2) * (m_q**2) / ( (1.0+m_q)**4 )    \
                                    + (   537.0*(m_q**6) - 3497.0*(m_q**5) - 18707.0*(m_q**4)      \
                                        - 29361.0*(m_q**3) - 18707.0*(m_q**2) - 3497.0*m_q + 537.0 \
                                      ) * (m_q/384.0) / ( (1.0+m_q)**8 )                           \
                                    - (1.0/16.0) * (S2[0]**2)    * (m_q**3) * ( 52.0*(m_q**2) + 12.0*m_q - 25.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/8.0)  * (S2[1]**2)    * (m_q**3) * (       m_q**2  - 17.0*m_q - 15.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S2[2]**2)    * (m_q**3) * ( 50.0*(m_q**2) + 38.0*m_q + 3.0  ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S1[0]**2)    * m_q      * ( 25.0*(m_q**2) - 12.0*m_q - 52.0 ) / ( (1.0+m_q)**6 )  \
                                    - (1.0/8.0)  * (S1[1]**2)    * m_q      * ( 15.0*(m_q**2) + 17.0*m_q -  1.0 ) / ( (1.0+m_q)**6 )  \
                                    + (1.0/16.0) * (S1[2]**2)    * m_q      * (  3.0*(m_q**2) + 38.0*m_q + 50.0 ) / ( (1.0+m_q)**6 )  \
                                    + (9.0/8.0)  * S1[0] * S2[0] * (m_q**3)                                       / ( (1.0+m_q)**6 )  \
                                    - (3.0/4.0)  * S1[1] * S2[1] * (m_q**2) * (  4.0*(m_q**2) +  9.0*m_q + 4.0)   / ( (1.0+m_q)**6 )  \
                                    + (3.0/8.0)  * S1[2] * S2[2] * (m_q**2) * ( 10.0*(m_q**2) + 21.0*m_q + 10.0 ) / ( (1.0+m_q)**6 )  \
                                  )

        dADM_dr_0PN  = mass * (   dADM_correction_0PN )
        dADM_dr_1PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN )
        dADM_dr_15PN = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN 
                              )
        dADM_dr_2PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                              )
        dADM_dr_25PN = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                                + dADM_correction_25PN                           \
                              )
        dADM_dr_3PN  = mass * (   dADM_correction_0PN  + dADM_correction_1PN      \
                                + dADM_correction_15PN + dADM_correction_2PN      \
                                + dADM_correction_25PN + dADM_correction_3PN      \
                              )
        
        return dADM_dr_0PN, dADM_dr_1PN, dADM_dr_15PN, dADM_dr_2PN, dADM_dr_25PN, dADM_dr_3PN
    
    ############################
    
    ADM_Mass_0PN,  ADM_Mass_1PN, ADM_Mass_15PN, ADM_Mass_2PN, ADM_Mass_25PN, ADM_Mass_3PN  = M_ADM(R)

    print(                                                           )
    print( " 0   阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_0PN  )
    print( " 1   阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_1PN  )
    print( " 1.5 阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_15PN )
    print( " 2   阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_2PN  )
    print( " 2.5 阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_25PN )
    print( " 3   阶后牛顿近似下的 ADM 质量 ADM Mass = ", ADM_Mass_3PN  )
    print(                                                           )

    Omega = Omega_3PN
    ADM_Mass_another_0PN,  ADM_Mass_another_1PN,  ADM_Mass_another_15PN, \
    ADM_Mass_another_2PN,  ADM_Mass_another_25PN, ADM_Mass_another_3PN   \
    = M_ADM_another(Omega)

    print(                                                                                         )
    print( " 0   阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_0PN  )
    print( " 1   阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_1PN  )
    print( " 1.5 阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_15PN )
    print( " 2   阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_2PN  )
    print( " 2.5 阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_25PN )
    print( " 3   阶后牛顿近似下的 ADM 质量（用 Omega 展开式计算） ADM Mass = ", ADM_Mass_another_3PN  )
    print(                                                                                         )

    ############################

    ## 利用设定的数值差分方法计算轨道能量 H 对 r 的微分
    ## 根据表达式 M_adm = M + H_circular               
    ##           ( 其中 M = M1 + M2 )

    def dH_dr(r):
        
        # 使用差分方法数值求导
        dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN \
                 = derivative.first_order_derivative_multivalue( M_ADM, r, 0.05, "7-points 6-orders" )

        # 使用 sympy 符号求导
        # 报错
        # dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN \
        #          = sympy.diff(M_ADM(r), r)
        
        return dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN
    
    ############################

    # 不知道为什么，这个差分的精度非常差
    # 可能是遭遇了舍入误差
    '''
    dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN = dH_dr(R)

    print(                                            )
    print( " 0   阶后牛顿近似下的 dH/dr = ", dHdr_0PN  )
    print( " 1   阶后牛顿近似下的 dH/dr = ", dHdr_1PN  )
    print( " 1.5 阶后牛顿近似下的 dH/dr = ", dHdr_15PN )
    print( " 2   阶后牛顿近似下的 dH/dr = ", dHdr_2PN  )
    print( " 2.5 阶后牛顿近似下的 dH/dr = ", dHdr_25PN )
    print(                                            )
    '''

    ############################

    ## 利用解析表达式直接计算轨道能量 H 对 r 的微分
    ## 根据表达式 M_adm = M + H_circular               
    ##           ( 其中 M = M1 + M2 )

    dHdr_0PN, dHdr_1PN, dHdr_15PN, dHdr_2PN, dHdr_25PN, dHdr_3PN = dADM_dr(R)

    print(                                            )
    print( " 0   阶后牛顿近似下的 dH/dr = ", dHdr_0PN  )
    print( " 1   阶后牛顿近似下的 dH/dr = ", dHdr_1PN  )
    print( " 1.5 阶后牛顿近似下的 dH/dr = ", dHdr_15PN )
    print( " 2   阶后牛顿近似下的 dH/dr = ", dHdr_2PN  )
    print( " 2.5 阶后牛顿近似下的 dH/dr = ", dHdr_25PN )
    print( " 3   阶后牛顿近似下的 dH/dr = ", dHdr_3PN  )
    print(                                            )

    ########################################

    ## 下面求出双星轨道的轨道间距随时间的变化率
    ## 根据文献  
    ## Antoni Ramos-Buades, Sascha Husa, and Geraint Pratten
    ## Simple procedures to reduce eccentricity of binary black hole simulations 
    ## arXiv:1810.00036[gr-qc]
    ## Phys. Rev.D 99, 02300 (2019) 

    ## 该文献关于 dE_GW/dt 的高阶结果有疑问
    ## 这里使用另一篇文献的结果
    ## Serguei Ossokine, Michael Boyle, Lawrence E. Kidder, Harald P. Pfeiffer, Mark A. Scheel, and Bela Szilagyi
    ## Comparing Post-Newtonian and Numerical-Relativity Precession Dynamics
    ## arXiv:1502.01747[gr-qc]
    ## Phys.Rev.D 92, 104028 (2015)

    ############################


    ## 取双星轨道运动的角速度为 3 阶后牛顿近似取值
    ## 各阶后牛顿近似取值前面已经算出

    Omega = Omega_3PN

    ## 算出引力波随时间变化率

    ## 根据文献 arXiv:1502.01747[gr-qc] 和 arXiv:1810.00036[gr-qc] 的定义，设置以下参数
    ## eta = m_eta
    m_delta     = ( m1 - m2 ) / ( m1 + m2 )
    
    spin_chi_a_vector = (S1 - S2) / 2.0
    spin_chi_s_vector = (S1 + S2) / 2.0
    
    spin_chi_a_square = numpy.dot(spin_chi_a_vector, spin_chi_a_vector)
    spin_chi_s_square = numpy.dot(spin_chi_s_vector, spin_chi_s_vector)

    ## 这里将向量 l 指向 z 轴正方向
    spin_chi_a_l  = ( S1[2] - S2[2]) / 2.0 
    spin_chi_s_l  = ( S1[2] + S2[2]) / 2.0

    Euler_gamma   = 0.5772156649

    mass    = m1 + m2
    Spin_l  = ( S1[2]*(m1**2) + S2[2]*(m2**2) ) / (mass**2)
    Sigma_l = ( m2*S2[2] - m1*S1[2] ) / mass

    ## 在 arXiv:1810.00036[gr-qc] 的公式中
    ## 前面需要手动加个负号
    dEGW_dt_0PN = - ( 32.0 / 5.0 ) * (m_eta**2) * ( Omega**(10.0/3.0) )

    dEGW_dt_correction_1PN  = - ( Omega**(2.0/3.0) ) * ( 35.0*m_eta/12.0 + 1247.0/336.0 )          \
                              + Omega * ( 4.0*math.pi - (5.0/4.0)*m_delta*Sigma_l - 4.0*Spin_l ) 

    dEGW_dt_correction_15PN = Omega**(4.0/3.0) \
                              * (   (65.0/18.0)    * (m_eta**2)   \
                                  + (9271.0/504.0) * m_eta        \
                                  - (44711.0/9072.0)              \
                                  -  (89.0/48.0) * m_delta * numpy.dot(spin_chi_a_vector, spin_chi_s_vector) \
                                  + (287.0/48.0) * m_delta * spin_chi_a_l * spin_chi_s_l                     \
                                  + (    287.0/96.0 -       12.0*m_eta )  * (spin_chi_a_l**2)                \
                                  + (    287.0/96.0 + (1.0/24.0)*m_eta )  * (spin_chi_s_l**2)                \
                                  + ( - (89.0/96.0) +        4.0*m_eta )  * spin_chi_a_square                \
                                  - (    89.0/96.0  + (7.0/24.0)*m_eta )  * spin_chi_s_square                \
                                )
    dEGW_dt_correction_2PN_part1 = Omega**(5.0/3.0) \
                                    * ( - math.pi           * ( (583.0/24.0)*m_eta + 8191.0/672.0 ) \
                                        + m_delta * Sigma_l * (   (43.0/4.0)*m_eta - (13.0/16.0)  ) \
                                        + Spin_l            * (  (272.0/9.0)*m_eta - 4.5          ) \
                                      ) 

    ## 以下是文献 arXiv:1502.01747[gr-qc] 的结果，与但文献 arXiv:1810.00036[gr-qc] 的结果并不相同                                                           
    dEGW_dt_correction_2PN_part2 = (Omega**2) \
                                    * ( -    (775.0/324.0) * (m_eta**3)                                \
                                        - (94403.0/3024.0) * (m_eta**2)                                \
                                        + ( - (134543.0/7776.0) + (41.0/48.0)*(math.pi**2) ) * m_eta   \
                                        + (6643739519.0/69854400.0)  + (16.0/3.0) * (math.pi**2)       \
                                        + (1712.0/105.0) * ( - Euler_gamma                             \
                                                             - 0.5*math.log(16.0*(Omega**(2.0/3.0)))   \
                                                           )                                           \
                                        - (31.0/6.0) * math.pi * m_delta * Sigma_l                     \
                                        -       16.0 * math.pi * Spin_l                                \
                                      )
    ## 以下各项是文献 arXiv:1810.00036[gr-qc] 的结果，但与文献 arXiv:1502.01747[gr-qc] 的结果并不相同
    dEGW_dt_correction_2PN_part2b = (Omega**2) \
                                    * ( - 4843497781.0/69854400.0                                  \
                                        - (775.0/324.0) * (m_eta**3)                               \
                                        - (94403.0/3024.0) * (m_eta**2)                            \
                                        + ( 8009293.0/54432.0 - (41.0/64.0)*(math.pi**2) ) * m_eta \
                                        + (287.0/192.0) * (math.pi**2)                             \
                                        + (1712.0/105.0) * ( - Euler_gamma                             \
                                                             + (35.0/107.0)*(math.pi**2)               \
                                                             - 0.5*math.log(16.0*(Omega**(2.0/3.0))) ) \
                                        - (31.0/6.0) * math.pi * m_delta * Spin_l                      \
                                        - 16.0 * math.pi * Spin_l                                      \
                                        + m_delta * spin_chi_a_l * spin_chi_s_l * (   611.0/252.0          \
                                                                                    - (809.0/18.0)*m_eta   \
                                                                                  )                        \
                                        + spin_chi_a_square * (   43.0*(m_eta**2)          \
                                                                - (8345.0/504.0)*m_eta     \
                                                                + 611.0/504.0              \
                                                              )                            \
                                        + spin_chi_s_square * (   (173.0/18.0)*(m_eta**2)  \
                                                                - (2393.0/72.0)*m_eta      \
                                                                + 611.0/504.0              \
                                                              )                            \
                                      )

    dEGW_dt_correction_25PN = Omega**(7.0/3.0)     \
                              * (   ( (1933585.0/3024.0)*(m_eta**2) + (214745.0/1728.0)*m_eta - (16258.0/504.0)   ) * math.pi           \
                                  + (    - (2810.0/27.0)*(m_eta**2) +    (6172.0/189.0)*m_eta + (476645.0/6784.0) ) * Spin_l            \
                                  + (    - (1501.0/36.0)*(m_eta**2) +    (1849.0/126.0)*m_eta + (9535.0/336.0)    ) * m_delta * Sigma_l \
                                )             

    dEGW_dt_correction_3PN  = Omega**(8.0/3.0) \
                              * (   ( - (7163.0/672.0) + (130583.0/2016.0) * m_eta ) * math.pi * m_delta * Sigma_l \
                                  + ( - (3485.0/96.0)  +    (13879.0/72.0) * m_eta ) * math.pi * Spin_l            \
                                )

    dEGW_dt_1PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN  )
    dEGW_dt_15PN = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN       \
                                       + dEGW_dt_correction_15PN 
                                 )
    dEGW_dt_2PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                 )
    dEGW_dt_25PN = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                       + dEGW_dt_correction_25PN       \
                                 )
    
    dEGW_dt_3PN  = dEGW_dt_0PN * ( 1.0 + dEGW_dt_correction_1PN        \
                                       + dEGW_dt_correction_15PN       \
                                       + dEGW_dt_correction_2PN_part1  \
                                       + dEGW_dt_correction_2PN_part2b \
                                       + dEGW_dt_correction_25PN       \
                                       + dEGW_dt_correction_3PN        \
                                 )
    
    print(                                   )
    print( " dEGW/dt 0pn   = ", dEGW_dt_0PN  )
    print( " dEGW/dt 1pn   = ", dEGW_dt_1PN  )
    print( " dEGW/dt 1.5pn = ", dEGW_dt_15PN )
    print( " dEGW/dt 2pn   = ", dEGW_dt_2PN  )
    print( " dEGW/dt 2.5pn = ", dEGW_dt_25PN )
    print( " dEGW/dt 3pn   = ", dEGW_dt_3PN  )
    print(                                   )

    ## 根据表达式 dr/dt = (dE_GW/dt) / (dH_circular/dr)
    ##           M_adm = M + H_circular               
    ##           ( 其中 M = M1 + M2 )

    drdt_0PN  = dEGW_dt_0PN  / dHdr_0PN 
    drdt_1PN  = dEGW_dt_1PN  / dHdr_1PN
    drdt_15PN = dEGW_dt_15PN / dHdr_15PN
    drdt_2PN  = dEGW_dt_2PN  / dHdr_2PN
    drdt_25PN = dEGW_dt_25PN / dHdr_25PN
    drdt_3PN  = dEGW_dt_3PN  / dHdr_3PN

    print(                                                               )
    print( " 0   阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_0PN  )
    print( " 1   阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_1PN  )
    print( " 1.5 阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_15PN )
    print( " 2   阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_2PN  )
    print( " 2.5 阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_25PN )
    print( " 3   阶后牛顿近似下双星旋转的径向速度 Vr = dr/dt = ", drdt_3PN  )
    print(                                                               )

    '''
    ## 旧的公式
    ## 根据文献  
    ## A. Gopakumar, Bala R. Iyer, and Sai Iyer 
    ## Second post-Newtonian gravitational radiation reaction for two-body systems: Nonspinning bodies 
    ## arXiv:gr-qc/9703075
    ## 精确度比较低
    drdt_0PN = - (64.0/5.0) * ( epsilon**3 ) * m_eta
    drdt_1PN = drdt_0PN * ( 1.0 - ((1751.0/336.0) + 7.0*m_eta/4.0) * epsilon )
    drdt_2PN = drdt_0PN * ( 1.0 - ((1751.0/336.0) + 7.0*m_eta/4.0) * epsilon                                \
                                + ( 303455.0/18144 + 40981.0*m_eta/2016.0 + m_eta**2.0/2.0 ) * (epsilon**2) \
                          )
    '''

    ## 设定轨道间距随时间的变化率为 1.5 阶精度（以后可以修改这个地方）

    drdt = drdt_3PN 
    print(" 双星轨道间距随时间的变化率 = ", drdt )
    
    ##########################

    ## 这里使用下列文献的公式
    ## James Healy, Carlos O. Lousto, Hiroyuki Nakano, and Yosef Zlochower
    ## Post-Newtonian Quasicircular Initial Orbits for Numerical Relativity
    ## arXiv:1702.00872[qr-qc]
    ## Class. Quant. Grav. 34, 145011 (2017)

    factor_0PN = (1.0+m_q)**2 / m_q

    factor_correction_1PN  = - 0.5 * epsilon * ( 7.0*(m_q**2) + 15.0*m_q + 7.0 ) / m_q
    factor_correction_15PN = 0.0
    factor_correction_2PN  = + (1.0/8.0) * (epsilon**2)                                                 \
                               * ( 47.0*(m_q**4) + 229.0*(m_q**3) + 363.0*(m_q**2) + 229.0*m_q + 47.0 ) \
                               / ( m_q * ( (1.0+m_q)**2 ) ) 
    factor_correction_25PN = + 0.25 * (epsilon**2.5)                                                    \
                               * (   ( 12.0*(m_q**2) + 11.0*m_q + 4.0  ) * S2[2] / (1.0 + m_q)          \
                                   + (  4.0*(m_q**2) + 11.0*m_q + 12.0 ) * S1[2] / ( (1.0 + m_q)*m_q )  \
                                 )
    factor_correction_3PN  = epsilon**3 \
                             * ( - (1.0/16.0) * ( (math.pi)**2 ) \
                                 - (1.0/48.0) * ( 363.0*(m_q**6) + 2608.0*(m_q**5) + 7324.0*(m_q**4)          \
                                                  + 10161.0*(m_q**3) + 7324.0*(m_q**2) + 2608.0*m_q + 363.0   \
                                                ) / ( m_q * ( (1.0+m_q)**4 ) )                                \
                                 + 0.25 * (S2[0]**2) * m_q * ( 18.0*(m_q**2) + 6.0*m_q +  5.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.75 * (S2[1]**2) * m_q * (  3.0*(m_q**2) +     m_q +  1.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.75 * (S2[2]**2) * m_q * (  3.0*(m_q**2) +     m_q +  1.0 ) / ( (1.0+m_q)**2 )           \
                                 + 0.25 * (S1[0]**2)       * (  5.0*(m_q**2) + 6.0*m_q + 18.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 - 0.75 * (S1[1]**2)       * (       m_q**2  +     m_q +  3.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 - 0.75 * (S1[2]**2)       * (       m_q**2  +     m_q +  3.0 ) / ( m_q * ( (1.0+m_q)**2 ) ) \
                                 +         S1[0] * S2[0]   * (  3.0*(m_q**2) -     m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.5  *  S1[1] * S2[1]   * (  3.0*(m_q**2) - 2.0*m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                                 - 0.5  *  S1[2] * S2[2]   * (  3.0*(m_q**2) - 2.0*m_q +  3.0 ) / ( (1.0+m_q)**2 )           \
                               )

    factor_1PN  = factor_0PN + factor_correction_1PN  
    factor_15PN = factor_0PN + factor_correction_1PN + factor_correction_15PN
    factor_2PN  = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN
    factor_25PN = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN + factor_correction_25PN
    factor_3PN  = factor_0PN + factor_correction_1PN + factor_correction_15PN  \
                             + factor_correction_2PN + factor_correction_25PN  \
                             + factor_correction_3PN

    Numerator = drdt  - (epsilon**3.5) * ( - 0.25 * S1[0] * S1[1] * m_q      * (      m_q + 6.0  ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S1[0] * S2[1] * (m_q**2) * (  6.0*m_q + 13.0 ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S2[0] * S1[1] * m_q      * ( 13.0*m_q + 6.0  ) / ( (1.0+m_q)**4 ) \
                                           - 0.25 * S2[0] * S2[1] * (m_q**2) * (  6.0*m_q + 1.0  ) / ( (1.0+m_q)**4 ) \
                                         )

    print(                                                           )
    print( " dr/dt = ", drdt                                         )
    print( " Numerator     in Pr calculation  3pn   = ", Numerator   )
    print( " devide factor in Pr calculation  0pn   = ", factor_0PN  )
    print( " devide factor in Pr calculation  1pn   = ", factor_1PN  )
    print( " devide factor in Pr calculation  1.5pn = ", factor_15PN )
    print( " devide factor in Pr calculation  2pn   = ", factor_2PN  )
    print( " devide factor in Pr calculation  2.5pn = ", factor_25PN )
    print( " devide factor in Pr calculation  3pn   = ", factor_3PN  )
    print(                                                           )

    Pr_0PN  = drdt_0PN  / factor_0PN
    Pr_1PN  = drdt_1PN  / factor_1PN
    Pr_15PN = drdt_15PN / factor_15PN
    Pr_2PN  = drdt_2PN  / factor_2PN
    Pr_25PN = drdt_25PN / factor_25PN
    Pr_3PN  = Numerator / factor_3PN

    print(                                                     )
    print( " 0   阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_0PN  )
    print( " 1   阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_1PN  )
    print( " 1.5 阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_15PN )
    print( " 2   阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_2PN  )
    print( " 2.5 阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_25PN )
    print( " 3   阶后牛顿近似下双星旋转的径向动量 Pr = ", Pr_3PN  )
    print(                                                     )

    ########################################

    ## 求出双星的动量，并写入文件
    
    ## 注意
    ## 这里为了和 AMSS-NCKU 的 TwoPuncture 程序输入相吻合
    ## 将较大质量黑洞设置在 y 轴正向，较小质量黑洞设置在 y 轴负向
    ## 动量满足 P1 = [ -|Pt|, -|Pr| ]
    ##         P2 = [ +|Pt|, +|Pr| ]
    ## 这样两个黑洞在 t=0 时位于 y 轴正向，做逆时针旋转

    momentum1[1] = - abs(Pr_3PN)
    momentum2[1] =   abs(Pr_3PN) 

    momentum1[0] = - abs(Pt_3PN) 
    momentum2[0] =   abs(Pt_3PN)  

    print(                                          )
    print(                                          )
    print(  " 双星径向动量大小 |Pr| = ", abs(Pr_3PN) )
    print(  " 双星切向动量大小 |Pt| = ", abs(Pt_3PN) )
    print(                                          )
    print(  " 双星在 t=0 时刻的动量为："              )
    print( f" P1 = {momentum1}  P2 = {momentum2}"   )
    print(                                          )

    print(                            )
    print( " 双星系统轨道动量设置完毕 " )
    print(                            )

    ##############################################################################################

    ## 将结果写入文件，供 Einstein Toolkit 和 AMSS-NCKU 程序使用

    # file1 = open( "BBH_parameter.output", "w" )
    file1 = open( os.path.join(input_data.File_directionary, "BBH_parameter.output"), "w")

    print(                                           file=file1 )
    print( " 双星系统轨道参数 ",                      file=file1 )
    print(                                           file=file1 )
    print( f" 双星质量：     M1 = {M1}  M2 = {M2} ",  file=file1 )
    print( " 无量纲质量比为： Q  = M1/M2 = ", m_q,    file=file1 )
    print( f" 双星无量纲自旋：S1 = {S1}  S2 = {S2} ", file=file1 )
    print(                                           file=file1 )
    print( " 双星在 t=0 时刻的坐标为：",               file=file1 )
    print( " X1 = ", position1[0],                   file=file1 ) 
    print( " Y1 = ", position1[1],                   file=file1 )
    print( " X2 = ", position2[0],                   file=file1 ) 
    print( " Y2 = ", position2[1],                   file=file1 )
    print(                                           file=file1 )
    print( " 双星在 t=0 时刻的动量为：",               file=file1 )
    print( " Pr  = ", Pr_3PN,                        file=file1 ) 
    print( " Pt  = ", Pt_3PN,                        file=file1 )
    print( " PX1 = - |Pt| = ", momentum1[0],         file=file1 )
    print( " PY1 = - |Pr| = ", momentum1[1],         file=file1 ) 
    print( " PX2 = + |Pt| = ", momentum2[0],         file=file1 )
    print( " PX2 = + |Pr| = ", momentum2[1],         file=file1 ) 
    print(                                           file=file1 )

    file1.close()
    
    return momentum1, momentum2

##############################################################################################


##############################################################################################

## 调用函数来计算轨道动量

## generate_BBH_orbit_parameters( M1, M2, S1, S2, D0, e0 )

##############################################################################################

