
##################################################################
##
## 该文件用于对 Puncture 数据进行初始化
## 小曲
## 2026/02/10 
## 2026/08/27 修改 
##
##################################################################


import numpy
import os 
import AMSS_NCKU_Input as input_data          ## 导入程序输入文件
import math

##################################################################



##################################################################

## 该函数用以导入 puncture 的信息

##################################################################

def generate_puncture_input_data( input_language ):


    if ( input_language == "Chinese" ):
        print(                                          )
        print( " 正在设定 Puncture 的位置、动量、角动量参数" )
        print(                                          )
    elif ( input_language == "English" ):
        print(                                                                )
        print( " Setting Puncture's position, momentum and angular momentum " )
        print(                                                                )
        

    #######################################################
    
    ## 如果黑洞坐标和动量的的设置方式为 TwoPuncture
    ## 则根据要求计算出初始的轨道坐标和轨道动量，并双黑洞将总质量 rescale 为 M=1

    if (input_data.puncture_data_set == "Automatically-BBH" ):

        mass_ratio_Q = input_data.parameter_BH[0,0] / input_data.parameter_BH[1,0]
    
        if ( mass_ratio_Q < 1.0 ):
            print( " 质量比设置错误，请重设！！！"   )
            print( " 将第一个黑洞设置为大质量！！！" )
        
        BBH_M1 = mass_ratio_Q / ( 1.0 + mass_ratio_Q )
        BBH_M2 = 1.0          / ( 1.0 + mass_ratio_Q )

        ## 导入双黑洞距离和偏心率
        distance_d0    = input_data.Distance_final
        ellipticity_e0 = input_data.e0
    
        ## 设置双黑洞的坐标
        ## 注意，这里自动调整，将较大质量黑洞放在 y 轴正向，将较小质量黑洞放在 y 轴负向
        ## TwoPuncture 程序输入需要有以下约定
        ## use Brugmann's convention
        ##  -----0-----> y
        ##   -      +     


        BBH_X1 = 0.0
        BBH_Y1 = distance_d0 * 1.0 / ( 1 + mass_ratio_Q )
        BBH_Z1 = 0.0

        BBH_X2 = 0.0
        BBH_Y2 = - distance_d0 * mass_ratio_Q / ( 1 + mass_ratio_Q )
        BBH_Z2 = 0.0
    
        position_BH    = numpy.zeros( (input_data.puncture_number,3) )
        position_BH[0] = [BBH_X1, BBH_Y1, BBH_Z1]
        position_BH[1] = [BBH_X2, BBH_Y2, BBH_Z2]
        
        ## 如果 puncture 数目超过两个，剩下的按照 input 文件中的数据直接读取
        for i in range(input_data.puncture_number):
            if (i>1):
                position_BH[i] = input_data.position_BH[i]
    
        ## 设置 AMSS-NCKU TwoPuncture 程序的自旋角动量输入
        ## 注意这里的角动量不是无量纲角动量，需要乘以质量平方才行
        ## 经过测试，这里需要乘的是比质量（即设置总质量为 1 ）
    
        ## angular_momentum_BH = input_data.angular_momentum_BH

        angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )  
        dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
        charge_Q_BH           = numpy.zeros( input_data.puncture_number ) 
    
        for i in range(input_data.puncture_number):
        
            ## 得到 puncture 的质量和电荷
            if i==0:
                dimensionless_mass_BH[i] = BBH_M1
                charge_Q_BH[i]           = BBH_M1 * input_data.parameter_BH[i,1]
            elif i==1:
                dimensionless_mass_BH[i] = BBH_M2
                charge_Q_BH[i]           = BBH_M2 * input_data.parameter_BH[i,1]
            else:
                dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
                charge_Q_BH[i]           = input_data.parameter_BH[i,0] * input_data.parameter_BH[i,1]
    
            ## 得到 puncture 的自旋角动量
            
            if ( input_data.Symmetry == "octant-symmetry" ):
                if i==0:
                    angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M1**2) * input_data.parameter_BH[i,2] ]
                elif i==1:
                    angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M2**2) * input_data.parameter_BH[i,2] ]
                else:
                    angular_momentum_BH[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ] 
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
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
                    angular_momentum_BH[i] = (BBH_M2**2) * input_data.dimensionless_spin_BH[i]
                else:
                   angular_momentum_BH[i] = (input_data.parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i]
                   
        ## 设置两黑洞的无量纲自旋

        BBH_S1 = angular_momentum_BH[0] / BBH_M1**2
        BBH_S2 = angular_momentum_BH[1] / BBH_M2**2
    
        ## 从参数文件导入动量
        ## momentum_BH  = input_data.momentum_BH
    
        ## 从生成双星轨道动量的文件中计算轨道动量
        import BBH_orbit_parameter 

        momentum_BH = numpy.zeros( (input_data.puncture_number,3) )

        ## 根据后牛顿近似算出初始的轨道动量
        momentum_BH[0], momentum_BH[1] = BBH_orbit_parameter.generate_BBH_orbit_parameters( input_language,               \
                                                                                            BBH_M1,      BBH_M2,          \
                                                                                            BBH_S1,      BBH_S2,          \
                                                                                            distance_d0, ellipticity_e0 ) 
        
        ## 如果 puncture 数目超过两个，剩下的按照 input 文件中的数据直接读取
        for i in range(input_data.puncture_number):
            if (i>1):
                momentum_BH[i] = input_data.momentum_BH[i]
                   
            
    #######################################################

    ## 如果黑洞坐标和动量的的设置方式为 Manually
    ## 则直接从参数文件中读入初始的轨道坐标和轨道动量
    ## 这里同样将双黑洞将总质量 rescale 为 M=1 （因为这个文件是为 TwoPuncture 生成输入文件）

    elif (input_data.puncture_data_set == "Manually" ):
    
        parameter_BH = input_data.parameter_BH
        position_BH  = input_data.position_BH
        momentum_BH  = input_data.momentum_BH
        
        ## 导入双黑洞距离和偏心率
        distance_d0    = math.sqrt(   (position_BH[0,0]-position_BH[1,0])**2   \
                                    + (position_BH[0,1]-position_BH[1,1])**2   \
                                    + (position_BH[0,2]-position_BH[1,2])**2 )
        ellipticity_e0 = input_data.e0
        
        ## 不再将双黑洞总质量进行归一化，直接读入参数文件中的值
        
        angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )   
        dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
        charge_Q_BH           = numpy.zeros( input_data.puncture_number )
            
        for i in range(input_data.puncture_number):
                
            dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
            charge_Q_BH[i]           = dimensionless_mass_BH[i] * input_data.parameter_BH[i,1]
                
            if ( input_data.Symmetry == "octant-symmetry" ):
                angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ] 
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
                angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ]
            elif ( input_data.Symmetry == "no-symmetry" ):
                angular_momentum_BH[i] = (parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i]


        BBH_M1 = dimensionless_mass_BH[0]
        BBH_M2 = dimensionless_mass_BH[1]
                
                
        ## 在以下的设定中，如果 Puncture 的数目为两个，则自动将总质量归一化为 M=1
        
        ##############################
        '''
        ## 如果 Puncture 的数目为两个，则自动将总质量归一化为 M=1
    
        if (input_data.puncture_number == 2 ):

            mass_ratio_Q = input_data.parameter_BH[0,0] / input_data.parameter_BH[1,0]
    
            if ( mass_ratio_Q < 1.0 ):
                print( " 质量比设置错误，请重设！！！" )
                print( " 将第一个黑洞设置为大质量！！！" )
        
            BBH_M1 = mass_ratio_Q / ( 1.0 + mass_ratio_Q )
            BBH_M2 = 1.0          / ( 1.0 + mass_ratio_Q )
    
            ## 设置 AMSS-NCKU TwoPuncture 程序的自旋角动量输入
            ## 注意这里的角动量不是无量纲角动量，需要乘以质量平方才行
            ## 经过测试，这里需要乘的是比质量（即设置总质量为 1 ）

            ## angular_momentum_BH = input_data.angular_momentum_BH

            angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )   
            dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
            charge_Q_BH           = numpy.zeros( input_data.puncture_number ) 

        
            for i in range(input_data.puncture_number):
        
                ## 得到 puncture 的质量和电荷
                if i==0:
                    dimensionless_mass_BH[i] = BBH_M1
                    charge_Q_BH[i]           = BBH_M1 * input_data.parameter_BH[i,1]
                elif i==1:
                    dimensionless_mass_BH[i] = BBH_M2
                    charge_Q_BH[i]           = BBH_M2 * input_data.parameter_BH[i,1]
                else:
                    dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
                    charge_Q_BH[i]           = dimensionless_mass_BH[i] * input_data.parameter_BH[i,1]
            
                ## 得到 puncture 的自旋角动量
    
                if ( input_data.Symmetry == "octant-symmetry" ):
                    if i==0:
                        angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M1**2) * parameter_BH[i,2] ]
                    elif i==1:
                        angular_momentum_BH[i] = [ 0.0, 0.0, (BBH_M2**2) * parameter_BH[i,2] ]
                    else:
                        angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ] 
                elif ( input_data.Symmetry == "equatorial-symmetry" ):
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
                        
        
        ##############################
        
        ## 如果 Puncture 的数目不为2，则直接读入输入文件的参数，不将总质量归一化
    
        else:
            
            angular_momentum_BH   = numpy.zeros( (input_data.puncture_number, 3) )   
            dimensionless_mass_BH = numpy.zeros( input_data.puncture_number )
            charge_Q_BH           = numpy.zeros( input_data.puncture_number )
            
            for i in range(input_data.puncture_number):
                
                dimensionless_mass_BH[i] = input_data.parameter_BH[i,0]
                charge_Q_BH[i]           = dimensionless_mass_BH[i] * input_data.parameter_BH[i,1]
                
                if ( input_data.Symmetry == "octant-symmetry" ):
                    angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ] 
                elif ( input_data.Symmetry == "equatorial-symmetry" ):
                    angular_momentum_BH[i] = [ 0.0, 0.0, (parameter_BH[i,0]**2) * parameter_BH[i,2] ]
                elif ( input_data.Symmetry == "no-symmetry" ):
                    angular_momentum_BH[i] = (parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i]
        '''            
                    
    #######################################################
    
    ## 其它情况，直接报错
    
    else: 
   
        if ( input_language == "Chinese" ):
            print(                                             )
            print( " 设定 Puncture 位置、动量、角动量参数的设定错误" )
            print(                                             )
        elif ( input_language == "English" ):
            print(                                                                                    )
            print( " Found Error in setting Puncture's position, momentum, and angular momentum !!! " )
            print(                                                                                    )
            
    #######################################################
    
    ## 设定 class 作为结构体，用来返回变量

    class Puncture_Input_Class:
        def __init__( self, BBH_M1,                   BBH_M2,            \
                            distance_d0,              ellipticity_e0,    \
                            dimensionless_mass_BH,    charge_Q_BH,       \
                            position_BH,              momentum_BH,       \
                            angular_momentum_BH ):
            self.BBH_M1                = BBH_M1
            self.BBH_M2                = BBH_M2
            self.distance_d0           = distance_d0
            self.ellipticity_e0        = ellipticity_e0
            self.charge_Q_BH           = charge_Q_BH
            self.position_BH           = position_BH
            self.momentum_BH           = momentum_BH
            self.angular_momentum_BH   = angular_momentum_BH
            self.dimensionless_mass_BH = dimensionless_mass_BH
            
    ## 创建 class 实例
    puncture_data = Puncture_Input_Class( BBH_M1,                BBH_M2,         \
                                          distance_d0,           ellipticity_e0, \
                                          dimensionless_mass_BH, charge_Q_BH,    \
                                          position_BH,           momentum_BH,    \
                                          angular_momentum_BH )
    
    ## 以下为测试，检查是否创建正确的 class 实例
    '''
    print( puncture_data.BBH_M1                )
    print( puncture_data.BBH_M2                )
    print( puncture_data.distance_d0           )
    print( puncture_data.ellipticity_e0        )
    print( puncture_data.dimensionless_mass_BH )
    print( puncture_data.charge_Q_BH           )
    print( puncture_data.position_BH           )
    print( puncture_data.momentum_BH           )
    print( puncture_data.angular_momentum_BH   )
    print( puncture_data.dimensionless_mass_BH )
    '''
    
    #######################################################
    
    ## 函数结束
    
    if ( input_language == "Chinese" ):
        print(                                              )
        print( " 完成 Puncture 的位置、动量、角动量参数的设定" )
        print(                                              )
    elif ( input_language == "English" ):
        print(                                                                      )
        print( " Puncture's position, momentum and angular momentum have been set " )
        print(                                                                      )
        
    print_puncture_information( input_language, puncture_data )
    
    return puncture_data
    
##################################################################




##################################################################

## 该函数用于输出黑洞 puncture 的信息

def print_puncture_information( input_language, puncture_data ):

    mass               = numpy.zeros( input_data.puncture_number )              ## 初始化每个黑洞的质量
    charge             = numpy.zeros( input_data.puncture_number )              ## 初始化每个黑洞的电荷
    angular_momentum_a = numpy.zeros( input_data.puncture_number )              ## 初始化每个黑洞的无量纲自旋
    position           = numpy.zeros( (input_data.puncture_number, 3) )         ## 初始化每个黑洞的位置
    momentum           = numpy.zeros( (input_data.puncture_number, 3) )         ## 初始化每个黑洞的的动量
    angular_momentum   = numpy.zeros( (input_data.puncture_number, 3) )         ## 初始化每个黑洞的角动量
    parameter          = numpy.zeros( (input_data.puncture_number, 3) )         ## 初始化每个黑洞的的参数

    print("------------------------------------------------------------------------------------------") 
    
    if ( input_language == "Chinese" ):
        print(                                       )   
        print( " 下列输出黑洞 puncture 的信息 "        )
        print(                                       ) 
    elif ( input_language == "English" ):
        print(                                       )
        print( " Printing the puncture information " )
        print(                                       )   

    
    ## 根据输入文件的设置得到每个黑洞的参数
    
    for i in range(input_data.puncture_number):

        parameter[i] = input_data.parameter_BH[i]
        
        position[i]  = puncture_data.position_BH[i] 
        momentum[i]  = puncture_data.momentum_BH[i]
        
        angular_momentum[i] = puncture_data.angular_momentum_BH[i]
        
        mass[i]   = puncture_data.dimensionless_mass_BH[i]
        charge[i] = puncture_data.charge_Q_BH[i]

        ## 根据输入文件设置每个黑洞的真实角动量
        '''
        if ( input_data.Symmetry == "octant-symmetry" ):
            angular_momentum[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            angular_momentum[i] = [ 0.0, 0.0, (input_data.parameter_BH[i,0]**2) * input_data.parameter_BH[i,2] ]
        elif ( input_data.Symmetry == "no-symmetry" ):
            angular_momentum[i] = (input_data.parameter_BH[i,0]**2) * input_data.dimensionless_spin_BH[i] 
        '''      
        
        ## 得到黑洞的无量纲自旋
        if ( input_data.Symmetry == "octant-symmetry" ):
            angular_momentum_a[i] = angular_momentum[i,2] / ( mass[i]**2 )
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            angular_momentum_a[i] = angular_momentum[i,2] / ( mass[i]**2 )
        elif ( input_data.Symmetry == "no-symmetry" ):
            angular_momentum_a[i] = ( angular_momentum[i,0]**2 + angular_momentum[i,1]**2 + angular_momentum[i,2]**2 )**(0.5) / ( mass[i]**2 )

        
        if ( input_language == "Chinese" ):
            print( f" 第 {i+1} 个黑洞的信息 " )  
        elif ( input_language == "English" ):
            print( f" The information for puncture {i+1} " ) 
            
        print( f" Mass({i+1}) = {mass[i]              :>10.6f},  Q({i+1})  = {charge[i]            :>10.6f},  a*({i+1}) = {angular_momentum_a[i]:>10.6f}" )
        print( f" X({i+1})    = {position[i,0]        :>10.6f},  Y({i+1})  = {position[i,1]        :>10.6f},  Z({i+1})  = {position[i,2]        :>10.6f}" )
        print( f" Px({i+1})   = {momentum[i,0]        :>10.6f},  Py({i+1}) = {momentum[i,1]        :>10.6f},  Pz({i+1}) = {momentum[i,2]        :>10.6f}" )
        print( f" Jx({i+1})   = {angular_momentum[i,0]:>10.6f},  Jy({i+1}) = {angular_momentum[i,1]:>10.6f},  Jz({i+1}) = {angular_momentum[i,2]:>10.6f}" )
        print(                                                                                                                                            )   

    print("------------------------------------------------------------------------------------------") 
    
    ## 请用户检查 Puncture 信息是否合理，如果合理，按回车继续
    
    if ( input_language == "Chinese" ):
        print(                                                                            )
        print( " 检查 Puncture 信息是否满足合理 "                                            )
        print( " （包含了质量、电荷、无量纲角动量 a*、位置坐标、动量、；角动量）"                   )
        print(                                                                            )
        print( " 如果 Puncture 的信息不合理，Ctrl+C 退出，在输入文件中调整 Puncture 的设定！！！ " )
        print( " 如果 Puncture 的信息合理，按回车继续！！！ "                                   )
        print(                                                                            )
    elif ( input_language == "English" ):
        print(                                                                                           )
        print( " Please check whether the puncture parameters are appropriate "                          )
        print( " (including mass, charge, dimensionless spin a*, position, momentum, angular momentum) " )
        print(                                                                                           )
        print( " If the puncture parameters are not appropriate, press Ctrl+C to abort the simulation. " )
        print( " Change the parameter setteing in the input file !!! "                                   )
        print(                                                                                           )
        print( " If the puncture parameters are appropriate, press Enter to continue !!!  "              )
        print(                                                                                           )
        
    ## 如果输入文件中设定了更远距离进行后牛顿演化，输出相关信息
    if ( input_data.puncture_data_set == "Automatically-BBH" and input_data.Allow_PN_Evaluation == "yes" ):
        if ( input_language == "Chinese" ):
            print(                                                                                                         )
            print( " 以上 Puncture 的初始动量仅为估算，后续还会用 PNOrbit 后牛顿程序从更远的区域进行演化，更新 Puncture 的初始动量信息 " )
            print( " 后续 PNOrbit 程序将从 Puncture 间距进行后牛顿演化 d(t0) = ", input_data.Distance_initial                   )
            print( " 后续 PNOrbit 程序将后牛顿演化到间距 d(tf) = ",              input_data.Distance_final                     )
            print(                                                                                                         )
        elif ( input_language == "English" ):
            print(                                                                                                                                          )
            print( "The above initial momentum of the Puncture is only an estimate; "                                                                       )
            print( "Later, the PNOrbit post-Newtonian code will evolve it from a larger separation and update the Puncture's initial momentum information." )
            print( "The PNOrbit code will perform a post-Newtonian evolution starting from the Puncture separation d(t0) = ", input_data.Distance_initial   )
            print( "The PNOrbit code will evolve post-Newtonianly until the separation d(tf) = ",                             input_data.Distance_final     )
            print(                                                                                                                                          )
        
    ## 按回车继续
    ## inputvalue = input()           
    ## print()
    
    return
    
##################################################################



##################################################################

## 该函数根据 PNOrbit 的运行结果更新 Puncture 数据
## 仅当 puncture_data_set == "Automatically-BBH" 且 Allow_PN_Evaluation == "yes" 时调用
## 解析 PNorbit.dat，用演化结束时刻（final）的双黑洞参数更新 position、momentum、distance

def update_puncture_data_by_PNOrbit( input_language, puncture_data ):

    import parse_PNOrbit_output

    ## 屏幕输出正在解析 PNOrbit 结果
    if ( input_language == "Chinese" ):
        print(                                            )
        print( " 正在解析 AMSS-NCKU 程序 PNOrbit 的运行结果 " )
        print(                                            )
    elif ( input_language == "English" ):
        print(                                                                      )
        print( " Parsing the running results of AMSS-NCKU executable file PNOrbit " )
        print(                                                                      )

    ## PNorbit.dat 位于 AMSS_NCKU_output 目录下
    output_directory    = os.path.join( input_data.File_directory, "AMSS_NCKU_output" )
    PNOrbit_data_file   = os.path.join( output_directory,          "PNorbit.dat"      )

    ## 解析结果，返回起始和终止时刻的 puncture 数据，以及终止时刻的双星间距
    position_BH_initial, momentum_BH_initial, position_BH_final, momentum_BH_final, distance_d0 \
        = parse_PNOrbit_output.parse_PNOrbit_output( PNOrbit_data_file, output_directory )

    ## 用演化结束时刻的双黑洞参数更新 puncture_data（位置和动量）
    puncture_data.position_BH = position_BH_final
    puncture_data.momentum_BH = momentum_BH_final

    ## 同步更新双星间距 distance_d0（PNOrbit 演化结束时刻的实际间距）
    puncture_data.distance_d0 = distance_d0

    ## 屏幕输出更新后的 puncture 信息
    if ( input_language == "Chinese" ):
        print(                                                    )
        print( " 已根据 PNOrbit 演化结束时刻的结果更新 Puncture 数据 " )
        print(                                                    )
        print( " 更新后的黑洞位置和动量 "                             )
        print( f" X1 = {position_BH_final[0,0]:>10.6f},  Y1 = {position_BH_final[0,1]:>10.6f}"  )
        print( f" X2 = {position_BH_final[1,0]:>10.6f},  Y2 = {position_BH_final[1,1]:>10.6f}"  )
        print( f" PX1 = {momentum_BH_final[0,0]:>10.6f}, PY1 = {momentum_BH_final[0,1]:>10.6f}" )
        print( f" PX2 = {momentum_BH_final[1,0]:>10.6f}, PY2 = {momentum_BH_final[1,1]:>10.6f}" )
        print(                                                                                  )
    elif ( input_language == "English" ):
        print(                                                                     )
        print( " The puncture data has been updated by the PNOrbit final results " )
        print(                                                                     )
        print( " The updated positions and momenta of the black holes "            )
        print( f" X1 = {position_BH_final[0,0]:>10.6f},  Y1 = {position_BH_final[0,1]:>10.6f}"  )
        print( f" X2 = {position_BH_final[1,0]:>10.6f},  Y2 = {position_BH_final[1,1]:>10.6f}"  )
        print( f" PX1 = {momentum_BH_final[0,0]:>10.6f}, PY1 = {momentum_BH_final[0,1]:>10.6f}" )
        print( f" PX2 = {momentum_BH_final[1,0]:>10.6f}, PY2 = {momentum_BH_final[1,1]:>10.6f}" )
        print(                                                                                  )

    return puncture_data

##################################################################
    
    
