
##################################################################
##
## 该文件生成 AMSS-NCKU 程序后牛顿轨道模块 PNOrbit 的输入文件 | This file generates the input file of the PNOrbit post-Newtonian orbit module of the AMSS-NCKU program
## 小曲 | XiaoQu
## 2024/11/27
## 2026/08/27 修改 | Modified on 2026/08/27
##
##################################################################


import numpy
import os 
import AMSS_NCKU_Input as input_data          ## 导入程序输入文件 | Import the program input file
import math

##################################################################

## 单位约定（与 PN_orbit.C 一致）：G = c = 1，总质量 M = m1 + m2 = 1。 | Unit convention (consistent with PN_orbit.C): G = c = 1, total mass M = m1 + m2 = 1.
## 在此几何单位制下，角动量与质量的平方同量纲（[S] = M^2）。 | In this geometric unit system, angular momentum has the same dimension as the square of mass ([S] = M^2).

##################################################################

## 该函数用于生成 PNOrbit（后牛顿轨道演化程序）的输入文件 | This function generates the input file of PNOrbit (the post-Newtonian orbit evolution code)
##
## 参数说明： | Parameter description:
##   puncture_data  : 由 puncture_initialize.py 构建的双黑洞初始参数对象，包含以下 |   puncture_data  : the binary-black-hole initial-parameter object built by puncture_initialize.py, containing the following
##   BBH_M1, BBH_M2         : 两个黑洞的质量（m1 + m2 = 1） |   BBH_M1, BBH_M2         : masses of the two black holes (m1 + m2 = 1)
##   position_BH[i]         : 黑洞 i 的初始位置（大质量黑洞在 +y 侧） |   position_BH[i]         : initial position of black hole i (the more massive black hole on the +y side)
##   momentum_BH[i]         : 黑洞 i 的初始动量 |   momentum_BH[i]         : initial momentum of black hole i
##   angular_momentum_BH[i] : 黑洞 i 的自旋角动量 —— 物理自旋 S_a = chi_a * m_a^2 |   angular_momentum_BH[i] : spin angular momentum of black hole i -- physical spin S_a = chi_a * m_a^2
##   distance_d0            : 双星初始坐标间距 |   distance_d0            : initial coordinate separation of the binary
##
##################################################################
    
def generate_AMSSNCKU_PNOrbit_input( input_language, puncture_data ): 

    ## 生成输入文件 AMSS-NCKU-PNOrbit.input | Generate the input file AMSS-NCKU-PNOrbit.input
    ## （AMSS_NCKU_Program.py 随后把它复制并重命名为 PN_Orbit.par） | (AMSS_NCKU_Program.py later copies it and renames it PN_Orbit.par)
    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU-PNOrbit.input"), "w") 

    ## 以下被注释的 ABE:: 键是旧版直接生成 TwoPuncture (ABE) 输入的写法， | The commented-out ABE:: keys below are the old way of directly generating the TwoPuncture (ABE) input,
    ## 现已由 PNOrbit 后牛顿演化流程取代：演化得到的轨道参数最终仍会 | now replaced by the PNOrbit post-Newtonian evolution workflow: the orbital parameters obtained from the evolution will still
    ## 经 generate_TwoPuncture_input.py 写入 ABE:: 键。整块保留备查。 | be written into the ABE:: keys via generate_TwoPuncture_input.py. The whole block is kept for reference.
    '''
    print( "#  -----0-----> y",                                         file=file1 )
    print( "#   -      +      use Brugmann's convention",               file=file1 )
    print( "ABE::mp        = -1.0",                                     file=file1 )   ## 这里要写成负数，方便程序自动求解裸质量 | This must be written as a negative number so that the program can automatically solve for the bare mass
    print( "ABE::mm        = -1.0",                                     file=file1 )
    print( "# b            =  D/2",                                     file=file1 )
    print( "ABE::b         = ", ( puncture_data.distance_d0 / 2.0 ),    file=file1 )
    print( "ABE::P_plusx   = ", puncture_data.momentum_BH[0,0],         file=file1 )
    print( "ABE::P_plusy   = ", puncture_data.momentum_BH[0,1],         file=file1 )
    print( "ABE::P_plusz   = ", puncture_data.momentum_BH[0,2],         file=file1 )
    print( "ABE::P_minusx  = ", puncture_data.momentum_BH[1,0],         file=file1 )
    print( "ABE::P_minusy  = ", puncture_data.momentum_BH[1,1],         file=file1 )
    print( "ABE::P_minusz  = ", puncture_data.momentum_BH[1,2],         file=file1 )
    print( "ABE::S_plusx   = ", puncture_data.angular_momentum_BH[0,0], file=file1 )
    print( "ABE::S_plusy   = ", puncture_data.angular_momentum_BH[0,1], file=file1 )
    print( "ABE::S_plusz   = ", puncture_data.angular_momentum_BH[0,2], file=file1 )
    print( "ABE::S_minusx  = ", puncture_data.angular_momentum_BH[1,0], file=file1 )
    print( "ABE::S_minusy  = ", puncture_data.angular_momentum_BH[1,1], file=file1 )
    print( "ABE::S_minusz  = ", puncture_data.angular_momentum_BH[1,2], file=file1 )
    print( "ABE::Mp        = ", puncture_data.BBH_M1,                   file=file1 )
    print( "ABE::Mm        = ", puncture_data.BBH_M2,                   file=file1 )
    print( "ABE::admtol    =  1.e-8",                                   file=file1 )
    print( "ABE::Newtontol =  5.e-12",                                  file=file1 )
    print( "ABE::nA        =  50",                                      file=file1 )
    print( "ABE::nB        =  50",                                      file=file1 )
    print( "ABE::nphi      =  26",                                      file=file1 )
    print( "ABE::Newtonmaxit =  50",                                    file=file1 )
    '''

    ## 对称质量比 nu = m1*m2/M^2 | Symmetric mass ratio nu = m1*m2/M^2
    symmetric_mass_ratio = puncture_data.BBH_M1 * puncture_data.BBH_M2 / ( (puncture_data.BBH_M1 + puncture_data.BBH_M2)**2 )

    print( "PN::mp = ",                            puncture_data.BBH_M1,                   file=file1 )  
    print( "PN::mm = ",                            puncture_data.BBH_M2,                   file=file1 )
    print( "PN::symmetric mass ratio =",           symmetric_mass_ratio,                   file=file1 )
  
    print( "PN::S_plusx   = ",                     puncture_data.angular_momentum_BH[0,0] / ( puncture_data.BBH_M1**2 ), file=file1 )
    print( "PN::S_plusy   = ",                     puncture_data.angular_momentum_BH[0,1] / ( puncture_data.BBH_M1**2 ), file=file1 )
    print( "PN::S_plusz   = ",                     puncture_data.angular_momentum_BH[0,2] / ( puncture_data.BBH_M1**2 ), file=file1 )
    print( "PN::S_minusx  = ",                     puncture_data.angular_momentum_BH[1,0] / ( puncture_data.BBH_M2**2 ), file=file1 )
    print( "PN::S_minusy  = ",                     puncture_data.angular_momentum_BH[1,1] / ( puncture_data.BBH_M2**2 ), file=file1 )
    print( "PN::S_minusz  = ",                     puncture_data.angular_momentum_BH[1,2] / ( puncture_data.BBH_M2**2 ), file=file1 )

    ## angular_momentum_BH 是"物理自旋"（有量纲量），不是无量纲自旋： | angular_momentum_BH is the "physical spin" (a dimensionful quantity), not the dimensionless spin:
    ## 而 PN_Orbit_Spin 版 C 程序的参数文件按文献惯例读取"无量纲自旋" chi_a | while the parameter file of the PN_Orbit_Spin version of the C code reads the "dimensionless spin" chi_a following the literature convention

    print( "PN::initial coordinate separation =" , input_data.Distance_initial, file=file1 )
    print( "PN::wanted r = ",                      input_data.Distance_final,   file=file1 )
    print( "PN::dT = 0.05",                                                     file=file1 )
    print( "PN::dumptime = 5.0",                                                file=file1 )
    print( "PN::totaltime = 10000000.0",                                        file=file1 )
    
    file1.close()

    return 
    
##################################################################
    
    

