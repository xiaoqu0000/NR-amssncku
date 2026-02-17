
##################################################################
##
## 该文件设置 AMSS-NCKU 程序的 TwoPuncture 的输入文件
## 小曲
## 2024/11/27
## 2026/02/10 修改
##
##################################################################


import numpy
import os 
import AMSS_NCKU_Input as input_data          ## 导入程序输入文件
import math

##################################################################




##################################################################

## 该函数用于生成 AMSS-NCKU-TwoPuncture 程序的输入文件

##################################################################
    
def generate_AMSSNCKU_TwoPuncture_input( input_language, puncture_data ): 

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU-TwoPuncture.input"), "w") 

    print( "#  -----0-----> y",                                         file=file1 )
    print( "#   -      +      use Brugmann's convention",               file=file1 )
    print( "ABE::mp        = -1.0",                                     file=file1 )   ## 这里要写成负数，方便程序自动求解裸质量
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
    
    file1.close()

    return 
    
##################################################################
    
    
