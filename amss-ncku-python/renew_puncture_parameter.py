
##################################################################
##
## 该文件根据 TwoPuncture 的输出结果来设置 Puncture 数据
## 小曲
## 2024/12/04
##
##################################################################

import AMSS_NCKU_Input as input_data
import numpy
import os

##################################################################



##################################################################

def read_TwoPuncture_Output(Output_File_directionary):

    dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
    bare_mass_BH          = numpy.zeros( input_data.puncture_number )        ## 初始化每个黑洞的质量
    position_BH           = numpy.zeros( (input_data.puncture_number, 3) )   ## 初始化每个黑洞的初始位置
    momentum_BH           = numpy.zeros( (input_data.puncture_number, 3) )   ## 初始化每个黑洞的动量
    angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )   ## 初始化每个黑洞的自旋角动量
    
    # 读取文件内容
    data = numpy.loadtxt( os.path.join(Output_File_directionary, "puncture_parameters_new.txt") )
    # 确保数据被解析为一维数组
    data = data.reshape(-1)
    ## print(" 读取到的 TwoPuncture 数据为 ")
    ## print(data)
    
    for i in range(input_data.puncture_number):
        
        ## 从 Two Puncture 的输出中读取前两个黑洞的参数
        ## 后面黑洞的参数从输入文件读取
        if i<2:
            bare_mass_BH[i]          = data[12*i]
            dimensionless_mass_BH[i] = data[12*i+1]
            position_BH[i]           = [ data[12*i+3], data[12*i+4],  data[12*i+5]  ]
            momentum_BH[i]           = [ data[12*i+6], data[12*i+7],  data[12*i+8]  ]
            angular_momentum_BH[i]   = [ data[12*i+9], data[12*i+10], data[12*i+11] ]
        else:
            dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
            bare_mass_BH[i]          = input_data.parameter_BH[i,0]
            position_BH[i]           = input_data.position_BH[i]
            momentum_BH[i]           = input_data.momentum_BH[i]
            ## 根据对称性读入角动量
            if ( input_data.Symmetry == "equatorial-symmetry" ):
                angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "no-symmetry" ):
                angular_momentum_BH[i] = (dimensionless_mass_BH[i]**2) * input_data.dimensionless_spin_BH[i]
    
    return bare_mass_BH, dimensionless_mass_BH, position_BH, momentum_BH, angular_momentum_BH
    
##################################################################


##################################################################

## 将以上格点信息追加写入到 AMSS-NCKU-TwoPuncture 程序的输入文件

def append_AMSSNCKU_BSSN_input(File_directionary, TwoPuncture_File_directionary): 

    charge_Q_BH = numpy.zeros( input_data.puncture_number )   ## 初始化每个黑洞的电荷

    ##  如果用 Ansorg-TwoPuncture 求解数值相对论初值，则从 TwoPuncture 的计算结果中读取裸质量、位置、角动量等参数
    if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
        bare_mass_BH, dimensionless_mass_BH, position_BH, momentum_BH, angular_momentum_BH = read_TwoPuncture_Output(TwoPuncture_File_directionary)
        # 设置每个黑洞电荷
        for i in range(input_data.puncture_number):
            charge_Q_BH[i] = dimensionless_mass_BH[i] * input_data.parameter_BH[i,1]
    
    ## 如果用其它方式求解数值相对论初值，则从输入文件直接读入参数    
    else:
        position_BH = input_data.position_BH
        momentum_BH = input_data.momentum_BH
        ## angular_momentum_BH = input_data.angular_momentum_BH

        angular_momentum_BH = numpy.zeros( (input_data.puncture_number, 3) )   ## 初始化每个黑洞的自旋角动量
        mass_BH             = numpy.zeros( input_data.puncture_number      )   ## 初始化每个黑洞的质量

        ## 设置每个黑洞的电荷和自旋角动量
        for i in range(input_data.puncture_number):

            if ( input_data.Symmetry == "octant-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                charge_Q_BH[i]         = mass_BH[i]* input_data.parameter_BH[i,1]
                angular_momentum_BH[i] = [ 0.0, 0.0, (mass_BH[i]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                charge_Q_BH[i]         = mass_BH[i]* input_data.parameter_BH[i,1]
                angular_momentum_BH[i] = [ 0.0, 0.0, (mass_BH[i]**2) * input_data.parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "no-symmetry" ):
                mass_BH[i]             = input_data.parameter_BH[i,0]
                angular_momentum_BH[i] = (mass_BH[i]**2) * input_data.dimensionless_spin_BH[i]
                charge_Q_BH[i]         = mass_BH[i]      * input_data.parameter_BH[i,1]

    file1 = open( os.path.join(input_data.File_directionary, "AMSS-NCKU.input"), "a")   ## "a" 表示追加输出

    ## 输出 BSSN 的相关设定
    
    print(                                                                           file=file1 )
    print( "BSSN::chitiny  = 1e-5",                                                  file=file1 ) 
    print( "BSSN::time refinement start from level = ", input_data.refinement_level, file=file1 )
    print( "BSSN::BH_num   =  ",                        input_data.puncture_number,  file=file1 )
    
    for i in range(input_data.puncture_number):
    
        if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
            print( f"BSSN::Mass[{i}]  = { bare_mass_BH[i] } ",      file=file1 )
        else:
            print( f"BSSN::Mass[{i}]  = { mass_BH[i] } ",           file=file1 )
            
        print( f"BSSN::Qchar[{i}] = { charge_Q_BH[i] } ",           file=file1 )
        print( f"BSSN::Porgx[{i}] = { position_BH[i,0] } ",         file=file1 )
        print( f"BSSN::Porgy[{i}] = { position_BH[i,1] } ",         file=file1 )
        print( f"BSSN::Porgz[{i}] = { position_BH[i,2] } ",         file=file1 )
        print( f"BSSN::Pmomx[{i}] = { momentum_BH[i,0] } ",         file=file1 )
        print( f"BSSN::Pmomy[{i}] = { momentum_BH[i,1] } ",         file=file1 )
        print( f"BSSN::Pmomz[{i}] = { momentum_BH[i,2] } ",         file=file1 )
        print( f"BSSN::Spinx[{i}] = { angular_momentum_BH[i,0] } ", file=file1 )
        print( f"BSSN::Spiny[{i}] = { angular_momentum_BH[i,1] } ", file=file1 )
        print( f"BSSN::Spinz[{i}] = { angular_momentum_BH[i,2] } ", file=file1 )
            
    print(                                                          file=file1 )
    
    file1.close()

    return
    
#################################################

