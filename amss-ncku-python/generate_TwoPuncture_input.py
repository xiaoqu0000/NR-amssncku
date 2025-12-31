
##################################################################
##
## 该文件设置 AMSS-NCKU 程序的 TwoPuncture 的输入文件
## 小曲
## 2024/11/27
## 2025/01/21 修改
##
##################################################################


import numpy
import os 
import AMSS_NCKU_Input as input_data          ## 导入程序输入文件
import math

##################################################################

## 导入双黑洞的坐标

## 如果黑洞坐标和动量的的设置方式为 TwoPuncture
## 则根据要求计算出初始的轨道坐标和轨道动量，并双黑洞将总质量 rescale 为 M=1

if (input_data.puncture_data_set == "Automatically-BBH" ):

    mass_ratio_Q = input_data.parameter_BH[0,0] / input_data.parameter_BH[1,0]
    
    if ( mass_ratio_Q < 1.0 ):
        print( " 质量比设置错误，请重设！！！" )
        print( " 将第一个黑洞设置为大质量！！！" )
        
    BBH_M1 = mass_ratio_Q / ( 1.0 + mass_ratio_Q )
    BBH_M2 = 1.0          / ( 1.0 + mass_ratio_Q )

    ## 导入双黑洞距离和偏心率
    distance = input_data.Distance
    e0       = input_data.e0
    
    ## 设置双黑洞的坐标
    ## 注意，这里自动调整，将较大质量黑洞放在 y 轴正向，将较小质量黑洞放在 y 轴负向
    ## TwoPuncture 程序输入需要有以下约定
    ## use Brugmann's convention
    ##  -----0-----> y
    ##   -      +     


    BBH_X1 = 0.0
    BBH_Y1 = distance * 1.0 / ( 1 + mass_ratio_Q )
    BBH_Z1 = 0.0

    BBH_X2 = 0.0
    BBH_Y2 = - distance * mass_ratio_Q / ( 1 + mass_ratio_Q )
    BBH_Z2 = 0.0
    
    position_BH    = numpy.zeros( (2,3) )
    position_BH[0] = [BBH_X1, BBH_Y1, BBH_Z1]
    position_BH[1] = [BBH_X2, BBH_Y2, BBH_Z2]
    
    ## 从参数文件导入动量
    ## momentum_BH  = input_data.momentum_BH
    
    ## 从生成双星轨道动量的文件中计算轨道动量
    import BBH_orbit_parameter 

    ## 带入 BBH_orbit_parameter 中设置的双黑洞无量纲角动量
    BBH_S1 = BBH_orbit_parameter.S1
    BBH_S2 = BBH_orbit_parameter.S2

    momentum_BH = numpy.zeros( (2,3) )

    ## 根据后牛顿近似算出初始的轨道动量
    momentum_BH[0], momentum_BH[1] = BBH_orbit_parameter.generate_BBH_orbit_parameters( BBH_M1, BBH_M2, BBH_S1, BBH_S2, distance, e0 ) 
    
    ## 设置 AMSS-NCKU TwoPuncture 程序的自旋角动量输入
    ## 注意这里的角动量不是无量纲角动量，需要乘以质量平方才行
    ## 经过测试，这里需要乘的是比质量（即设置总质量为 1 ）
    
    ## angular_momentum_BH = input_data.angular_momentum_BH

    angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )  
    
    for i in range(input_data.puncture_number):
    
        if ( input_data.Symmetry == "equatorial-symmetry" ):
            if i==0:
                angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M1**2) * input_data.parameter_BH[i,2] ]
            elif i==1:
                angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M2**2) * input_data.parameter_BH[i,2] ]
            else:
                angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
                
        elif ( input_data.Symmetry == "no-symmetry" ):
        
            if i==0:
                angular_momentum_BH[i] = (BBH_M1**2) * input_data.dimensionless_spin_BH[i]
            elif i==1:
                angular_momentum_BH[i] = (BBH_M1**2) * input_data.dimensionless_spin_BH[i]
            else:
                angular_momentum_BH[i] = (input_data.parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i]
            
    #######################################################

## 如果黑洞坐标和动量的的设置方式为 Manually
## 则直接从参数文件中读入初始的轨道坐标和轨道动量
## 这里同样将双黑洞将总质量 rescale 为 M=1 （因为这个文件是为 TwoPuncture 生成输入文件）

elif (input_data.puncture_data_set == "Manually" ):

    mass_ratio_Q = input_data.parameter_BH[0,0] / input_data.parameter_BH[1,0]
    
    if ( mass_ratio_Q < 1.0 ):
        print( " 质量比设置错误，请重设！！！" )
        print( " 将第一个黑洞设置为大质量！！！" )
        
    BBH_M1 = mass_ratio_Q / ( 1.0 + mass_ratio_Q )
    BBH_M2 = 1.0          / ( 1.0 + mass_ratio_Q )
    
    parameter_BH = input_data.parameter_BH
    position_BH  = input_data.position_BH
    momentum_BH  = input_data.momentum_BH
    
    ## 导入双黑洞距离和偏心率
    distance = math.sqrt( (position_BH[0,0]-position_BH[1,0])**2 + (position_BH[0,1]-position_BH[1,1])**2 + (position_BH[0,2]-position_BH[1,2])**2 )
    e0       = input_data.e0
    
    ## 设置 AMSS-NCKU TwoPuncture 程序的自旋角动量输入
    ## 注意这里的角动量不是无量纲角动量，需要乘以质量平方才行
    ## 经过测试，这里需要乘的是比质量（即设置总质量为 1 ）

    ## angular_momentum_BH = input_data.angular_momentum_BH

    angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )   

        
    for i in range(input_data.puncture_number):
    
        if ( input_data.Symmetry == "equatorial-symmetry" ):
            if i==0:
                angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M1**2) * parameter_BH[i,2] ]
            elif i==1:
                angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M2**2) * parameter_BH[i,2] ]
            else:
                angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ]
                
        elif ( input_data.Symmetry == "no-symmetry" ):
            if i==0:
                angular_momentum_BH[i] = (BBH_M1**2) * input_data.dimensionless_spin_BH[i]
            elif i==1:
                angular_momentum_BH[i] = (BBH_M2**2) * input_data.dimensionless_spin_BH[i]
            else:
                angular_momentum_BH[i] = (parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i]


##################################################################

## 将以上格点信息写入 AMSS-NCKU-TwoPuncture 程序的输入文件
    
def generate_AMSSNCKU_TwoPuncture_input(): 

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU-TwoPuncture.input"), "w") 

    print( "#  -----0-----> y",                           file=file1 )
    print( "#   -      +      use Brugmann's convention", file=file1 )
    print( "ABE::mp        = -1.0",                       file=file1 )   ## 这里要写成负数，方便程序自动求解裸质量
    print( "ABE::mm        = -1.0",                       file=file1 )
    print( "# b            =  D/2",                       file=file1 )
    print( "ABE::b         = ", ( distance / 2.0 ),       file=file1 )
    print( "ABE::P_plusx   = ", momentum_BH[0,0],         file=file1 )
    print( "ABE::P_plusy   = ", momentum_BH[0,1],         file=file1 )
    print( "ABE::P_plusz   = ", momentum_BH[0,2],         file=file1 )
    print( "ABE::P_minusx  = ", momentum_BH[1,0],         file=file1 )
    print( "ABE::P_minusy  = ", momentum_BH[1,1],         file=file1 )
    print( "ABE::P_minusz  = ", momentum_BH[1,2],         file=file1 )
    print( "ABE::S_plusx   = ", angular_momentum_BH[0,0], file=file1 )
    print( "ABE::S_plusy   = ", angular_momentum_BH[0,1], file=file1 )
    print( "ABE::S_plusz   = ", angular_momentum_BH[0,2], file=file1 )
    print( "ABE::S_minusx  = ", angular_momentum_BH[1,0], file=file1 )
    print( "ABE::S_minusy  = ", angular_momentum_BH[1,1], file=file1 )
    print( "ABE::S_minusz  = ", angular_momentum_BH[1,2], file=file1 )
    print( "ABE::Mp        = ", BBH_M1,                   file=file1 )
    print( "ABE::Mm        = ", BBH_M2,                   file=file1 )
    print( "ABE::admtol    =  1.e-8",                     file=file1 )
    print( "ABE::Newtontol =  5.e-12",                    file=file1 )
    print( "ABE::nA        =  50",                        file=file1 )
    print( "ABE::nB        =  50",                        file=file1 )
    print( "ABE::nphi      =  26",                        file=file1 )
    print( "ABE::Newtonmaxit =  50",                      file=file1 )
    
    file1.close()

    return file1
    
##################################################################
    
    
