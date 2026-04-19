
##################################################################
##
## 该文件定义程序的一些屏幕输出
## 小曲
## 2024/03/22
## 2026/02/27 修改
##
##################################################################

import AMSS_NCKU_Input as input_data
import numpy 
import os
import math

##################################################################


##################################################################

def setup( input_language ):

    devide_factor = input_data.devide_factor

    static_grid_level = input_data.static_grid_level
    moving_grid_level = input_data.moving_grid_level
    total_grid_level  = input_data.grid_level

    static_grid_number = input_data.static_grid_number
    moving_grid_number = input_data.moving_grid_number

    if ( input_data.Symmetry=="octant-symmetry" ):
        maximal_domain_size_static_x = numpy.array( [ 0.0, input_data.largest_box_xyz_max[0] ] )
        maximal_domain_size_static_y = numpy.array( [ 0.0, input_data.largest_box_xyz_max[1] ] )
        maximal_domain_size_static_z = numpy.array( [ 0.0, input_data.largest_box_xyz_max[2] ] )
    elif( input_data.Symmetry=="equatorial-symmetry" ):
        maximal_domain_size_static_x = numpy.array( [ input_data.largest_box_xyz_min[0], input_data.largest_box_xyz_max[0] ] )
        maximal_domain_size_static_y = numpy.array( [ input_data.largest_box_xyz_min[1], input_data.largest_box_xyz_max[1] ] )
        maximal_domain_size_static_z = numpy.array( [ 0.0, input_data.largest_box_xyz_max[2] ] )
    elif( input_data.Symmetry=="no-symmetry" ):
        maximal_domain_size_static_x = numpy.array( [ input_data.largest_box_xyz_min[0], input_data.largest_box_xyz_max[0] ] )
        maximal_domain_size_static_y = numpy.array( [ input_data.largest_box_xyz_min[1], input_data.largest_box_xyz_max[1] ] )
        maximal_domain_size_static_z = numpy.array( [ input_data.largest_box_xyz_min[2], input_data.largest_box_xyz_max[2] ] )
    else:
        if ( input_language == "Chinese" ):
            print(                       )   
            print( " 对称性设置错误！！ " )
            print(                       )
        elif ( input_language == "English" ):
            print(                                   )   
            print( " Symmetry Setting is Wrong !!! " )
            print(                                   )
    
    ## 固定网格为长方体网格
    minimal_domain_size_static_x =   maximal_domain_size_static_x / ( (devide_factor)**(static_grid_level-1) )
    minimal_domain_size_static_y =   maximal_domain_size_static_y / ( (devide_factor)**(static_grid_level-1) )
    minimal_domain_size_static_z =   maximal_domain_size_static_z / ( (devide_factor)**(static_grid_level-1) )
    
    ## 强行设定移动网格为立方体网格
    ## 这里要转化为 numpy 数组，否则报错
    maximal_domain_size_moving_x = numpy.array( [ - 0.5 * (minimal_domain_size_static_x[1] - minimal_domain_size_static_x[0]) / devide_factor,     \
                                                  + 0.5 * (minimal_domain_size_static_x[1] - minimal_domain_size_static_x[0]) / devide_factor  ] ) \
                                   * ( moving_grid_number / static_grid_number )   
    maximal_domain_size_moving = maximal_domain_size_moving_x
    minimal_domain_size_moving = maximal_domain_size_moving / ( (devide_factor)**(input_data.moving_grid_level-1) )
    ## 原来的设定
    ## maximal_domain_size_moving = ( minimal_domain_size_static_x / devide_factor ) * ( moving_grid_number / static_grid_number )

    ## 设定网格分辨率
    maximal_resolution_static = (input_data.largest_box_xyz_max[0] - input_data.largest_box_xyz_min[0]) / static_grid_number
    minimal_resolution_static = maximal_resolution_static / ( (devide_factor)**(static_grid_level-1) )
    maximal_resolution_moving = minimal_resolution_static / devide_factor
    minimal_resolution_moving = maximal_resolution_moving / ( (devide_factor)**(moving_grid_level-1) )

    ## 设定时间步长
    TimeStep = input_data.Courant_Factor * maximal_resolution_static / ( (devide_factor)**(input_data.refinement_level) )

    ## 设定 shell-patch 网格
    shell_grid_number              = input_data.shell_grid_number
    minimal_domain_size_shellpatch = input_data.largest_box_xyz_max[0]
    maximal_domain_size_shellpatch = input_data.largest_box_xyz_max[0] + maximal_resolution_static * shell_grid_number[2]
    shellpatch_resolution_R        = maximal_resolution_static
    shellpatch_resolution_theta    = 0.5 * math.pi / shell_grid_number[1]
    shellpatch_resolution_phi      = 0.5 * math.pi / shell_grid_number[0]
    
    ## 设定 class 作为结构体，用来返回变量
    
    class Setup_Data_Class:
        def __init__( self, minimal_domain_size_static_x,   minimal_domain_size_static_y,   minimal_domain_size_static_z, \
                            maximal_domain_size_static_x,   maximal_domain_size_static_y,   maximal_domain_size_static_z, \
                            minimal_domain_size_moving,     maximal_domain_size_moving,                                   \
                            minimal_resolution_static,      maximal_resolution_static,                                    \
                            minimal_resolution_moving,      maximal_resolution_moving,                                    \
                            minimal_domain_size_shellpatch, maximal_domain_size_shellpatch,                               \
                            shellpatch_resolution_R,        shellpatch_resolution_phi,      shellpatch_resolution_theta,  \
                            TimeStep ):                 
            self.minimal_domain_size_static_x   = minimal_domain_size_static_x
            self.minimal_domain_size_static_y   = minimal_domain_size_static_y
            self.minimal_domain_size_static_z   = minimal_domain_size_static_z
            self.maximal_domain_size_static_x   = maximal_domain_size_static_x
            self.maximal_domain_size_static_y   = maximal_domain_size_static_y
            self.maximal_domain_size_static_z   = maximal_domain_size_static_z
            self.minimal_domain_size_moving     = minimal_domain_size_moving 
            self.maximal_domain_size_moving     = maximal_domain_size_moving 
            self.minimal_resolution_static      = minimal_resolution_static   
            self.maximal_resolution_static      = maximal_resolution_static 
            self.minimal_resolution_moving      = minimal_resolution_moving   
            self.maximal_resolution_moving      = maximal_resolution_moving  
            self.minimal_domain_size_shellpatch = minimal_domain_size_shellpatch
            self.maximal_domain_size_shellpatch = maximal_domain_size_shellpatch     
            self.shellpatch_resolution_R        = shellpatch_resolution_R     
            self.shellpatch_resolution_phi      = shellpatch_resolution_phi   
            self.shellpatch_resolution_theta    = shellpatch_resolution_theta
            self.TimeStep                       = TimeStep
        
    ## 创建 class 实例
    setup_data = Setup_Data_Class( minimal_domain_size_static_x,   minimal_domain_size_static_y,   minimal_domain_size_static_z, \
                                   maximal_domain_size_static_x,   maximal_domain_size_static_y,   maximal_domain_size_static_z, \
                                   minimal_domain_size_moving,     maximal_domain_size_moving,                                   \
                                   minimal_resolution_static,      maximal_resolution_static,                                    \
                                   minimal_resolution_moving,      maximal_resolution_moving,                                    \
                                   minimal_domain_size_shellpatch, maximal_domain_size_shellpatch,                               \
                                   shellpatch_resolution_R,        shellpatch_resolution_phi,      shellpatch_resolution_theta,  \
                                   TimeStep )
                            
    
    return setup_data

##################################################################



##################################################################

## 这个函数用来输出整个程序的基本输入

def print_input_data_Chinese( setup_data, File_directory ):
    
    ## 屏幕输出
    ## 输出中文信息
    
    print( "------------------------------------------------------------------------------------------" ) 
    print(                                                            )
    print( " 下面输出本程序的输入参数 "                                 )
    print(                                                            )
    print( " AMSS-NCKU 程序调用的进程数目 = ", input_data.MPI_processes ) 
    print(                                                            )
    print( " 求解的方程/体系 = ", input_data.Equation_Class            )
    print( " 初值的设定      = ", input_data.Initial_Data_Method       )
    print(                                                            )
    print( " 起始演化时间 = ", input_data.Start_Evolution_Time         ) 
    print( " 最终演化时间 = ", input_data.Final_Evolution_Time         ) 
    print( " 最大迭代次数 = ", input_data.Evolution_Step_Number       )
    print( " Courant因子  = ", input_data.Courant_Factor              )
    print( " 演化耗散因子 = ", input_data.Dissipation                  )
    print( " 体系的对称性 = ", input_data.Symmetry                     )
    print( " 时间演化方法 = ", input_data.Time_Evolution_Method        ) 
    print( " 有限差分方法 = ", input_data.Finite_Diffenence_Method     )
    print(                                                            )
    print( " 固定网格的类型    = ", input_data.static_grid_type        )
    print( " 可移动网格的类型  = ", input_data.moving_grid_type        )
    print(                                                            )
    print( " 固定网格的层数    = ", input_data.static_grid_level       )      
    print( " 可移动网格的层数  = ", input_data.moving_grid_level       )      
    print( " 全体网格的层数    = ", input_data.grid_level              )
    print(                                                           )
    print( " 每层固定网格的格点数目   = ", input_data.static_grid_number )
    print( " 每层可移动网格的格点数目 = ", input_data.moving_grid_number )
    print(                                                              )
    print( " 最外层固定网格的范围 X 方向 = ", setup_data.maximal_domain_size_static_x )
    print( " 最外层固定网格的范围 Y 方向 = ", setup_data.maximal_domain_size_static_y )
    print( " 最外层固定网格的范围 Z 方向 = ", setup_data.maximal_domain_size_static_z  )
    print( " 最内层固定网格的范围 X 方向 = ", setup_data.minimal_domain_size_static_x  )
    print( " 最内层固定网格的范围 Y 方向 = ", setup_data.minimal_domain_size_static_y  )
    print( " 最内层固定网格的范围 Z 方向 = ", setup_data.minimal_domain_size_static_z  )
    
    if ( input_data.moving_grid_level > 0):
        print( " 最外层可移动网格的范围   = ", setup_data.maximal_domain_size_moving  )
        print( " 最内层可移动网格的范围   = ", setup_data.minimal_domain_size_moving  )

    print(                                                                       )
    print( " 最外层固定网格的分辨率   = ", setup_data.maximal_resolution_static   )
    print( " 最内层固定网格的分辨率   = ", setup_data.minimal_resolution_static   )
    
    if ( input_data.moving_grid_level > 0):
        print( " 最外层可移动网格的分辨率 = ", setup_data.maximal_resolution_moving   )
        print( " 最内层可移动网格的分辨率 = ", setup_data.minimal_resolution_moving   )
    
    print(                                                               )
    print( " 时间细化从第 ", input_data.refinement_level, " 层网格开始 "   )
    print( " 程序计算中的最粗时间步长 = ", setup_data.TimeStep             )
    print(                                                               )
    
    ## 如果网格类型设置为 Shell-Patch，输出 Shell-Patch 的网格信息
    if input_data.basic_grid_set == "Shell-Patch":
        print( " 程序计算中使用了 Shell-Patch 类型网格 "                         )
        print( " Shell-Patch 网格的格点数目 = ", input_data.shell_grid_number   )
        print( " Shell-Patch 网格的最小半径 = ", setup_data.minimal_domain_size_shellpatch )
        print( " Shell-Patch 网格的最大半径 = ", setup_data.maximal_domain_size_shellpatch )
        print( " Shell-Patch 网格径向分辨率 = ", setup_data.shellpatch_resolution_R        )
        print( " Shell-Patch 网格角向分辨率 = ", setup_data.shellpatch_resolution_phi,  \
                                                setup_data.shellpatch_resolution_theta    )
        print(                                                                 )
    elif input_data.basic_grid_set == "Patch":   
        print( " 程序计算中仅使用 Patch 类型网格，没用使用 Shell-Patch 类型网格 "  )
        print(                                                                 )
    else:
        print( " 网格类型设置错误！"                                             )
        print(                                                                  )

    print( "------------------------------------------------------------------------------------------" ) 
    
    ## 文件输出
    
    filepath = os.path.join( File_directory, "AMSS_NCKU_resolution" )
    file0    = open(filepath, 'w')
    
    print(                                                              file=file0 )
    print( " 下面输出本程序的输入参数 ",                                    file=file0 )
    print(                                                              file=file0 )
    print( " AMSS-NCKU 程序调用的进程数目 = ", input_data.MPI_processes, file=file0 ) 
    print(                                                             file=file0 )
    print( " 起始演化时间 = ", input_data.Start_Evolution_Time,         file=file0 ) 
    print( " 最终演化时间 = ", input_data.Final_Evolution_Time,         file=file0 ) 
    print( " Courant因子  = ", input_data.Courant_Factor,              file=file0 )
    print( " 演化耗散因子 = ", input_data.Dissipation,                  file=file0 )
    print( " 体系的对称性 = ", input_data.Symmetry,                     file=file0 )
    print( " 时间演化方法 = ", input_data.Time_Evolution_Method,        file=file0 ) 
    print( " 有限差分方法 = ", input_data.Finite_Diffenence_Method,     file=file0 )
    print(                                                             file=file0 )
    print( " 固定网格的类型    = ", input_data.static_grid_type,         file=file0 )
    print( " 可移动网格的类型  = ", input_data.moving_grid_type,         file=file0 )
    print(                                                             file=file0 )
    print( " 固定网格的层数    = ", input_data.static_grid_level,       file=file0 )      
    print( " 可移动网格的层数  = ", input_data.moving_grid_level,       file=file0 )      
    print( " 全体网格的层数    = ", input_data.grid_level,              file=file0 )
    print(                                                             file=file0 )
    print( " 每层固定网格的格点数目   = ", input_data.static_grid_number, file=file0 )
    print( " 每层可移动网格的格点数目 = ", input_data.moving_grid_number, file=file0 )
    print(                                                              file=file0 )
    print( " 最外层固定网格的范围 X 方向 = ", setup_data.maximal_domain_size_static_x, file=file0 )
    print( " 最外层固定网格的范围 Y 方向 = ", setup_data.maximal_domain_size_static_y, file=file0)
    print( " 最外层固定网格的范围 Z 方向 = ", setup_data.maximal_domain_size_static_z, file=file0  )
    print( " 最内层固定网格的范围 X 方向 = ", setup_data.minimal_domain_size_static_x, file=file0 )
    print( " 最内层固定网格的范围 Y 方向 = ", setup_data.minimal_domain_size_static_y, file=file0  )
    print( " 最内层固定网格的范围 Z 方向 = ", setup_data.minimal_domain_size_static_z, file=file0  )
    print(                                                                           file=file0 )
    print( " 最外层可移动网格的范围   = ", setup_data.maximal_domain_size_moving,     file=file0 )
    print( " 最内层可移动网格的范围   = ", setup_data.minimal_domain_size_moving,     file=file0 )
    print(                                                                           file=file0 )
    print( " 最外层固定网格的分辨率   = ", setup_data.maximal_resolution_static,      file=file0 )
    print( " 最内层固定网格的分辨率   = ", setup_data.minimal_resolution_static,      file=file0 )
    print( " 最外层可移动网格的分辨率 = ", setup_data.maximal_resolution_moving,      file=file0 )
    print( " 最内层可移动网格的分辨率 = ", setup_data.minimal_resolution_moving,      file=file0 )
    print(                                                                           file=file0 )
    print( " 时间细化从第 ", input_data.refinement_level, " 层网格开始 ",             file=file0 )
    print( " 程序计算中的最粗时间步长 = ", setup_data.TimeStep,                       file=file0 )
    print(                                                                           file=file0 )

    ## 如果网格类型设置为 Shell-Patch，输出 Shell-Patch 的网格信息
    if input_data.basic_grid_set == "Shell-Patch":
        print( " 程序计算中使用了 Shell-Patch 类型网格 ",                                    file=file0 )
        print( " Shell-Patch 网格的格点数目 = ", input_data.shell_grid_number,              file=file0 )
        print( " Shell-Patch 网格的最小半径 = ", setup_data.minimal_domain_size_shellpatch, file=file0 )
        print( " Shell-Patch 网格的最大半径 = ", setup_data.maximal_domain_size_shellpatch, file=file0 )
        print( " Shell-Patch 网格径向分辨率 = ", setup_data.shellpatch_resolution_R,        file=file0 )
        print( " Shell-Patch 网格角向分辨率 = ", setup_data.shellpatch_resolution_phi,  \
                                                setup_data.shellpatch_resolution_theta,    file=file0 )
        print(                                                                             file=file0 )
    elif input_data.basic_grid_set == "Patch":   
        print( " 程序计算中仅使用 Patch 类型网格，没用使用 Shell-Patch 类型网格 ",  file=file0 )
        print(                                                                  file=file0 )
    else:
        print( " 网格类型设置错误！",                                             file=file0 )
        print(                                                                  file=file0 )
        
    return
    
##################################################################    

## 这是同一份函数，用英语输出    
    
def print_input_data_English( setup_data, File_directory ):
    
    ## 屏幕输出  
    ## 下面输出英文信息
    
    print( "------------------------------------------------------------------------------------------" ) 
    print(                                                                                           )
    print(                                                                                           )
    print( " Printing the basic parameter and setting in the AMSS-NCKU simulation "                  )
    print(                                                                                           )
    print( " The number of MPI processes in the AMSS-NCKU simulation = ", input_data.MPI_processes   ) 
    print(                                                                                           )
    print( " The form of computational equation  = ",            input_data.Equation_Class           )
    print( " The initial data in this simulation = ",            input_data.Initial_Data_Method      )
    print(                                                                                           )
    print( " Starting evolution time   = ",                      input_data.Start_Evolution_Time     ) 
    print( " Final evolution time      = ",                      input_data.Final_Evolution_Time     ) 
    print( " Maximal iteration number  = ",                      input_data.Evolution_Step_Number    )
    print( " Courant factor            = ",                      input_data.Courant_Factor           )
    print( " Strength of dissipation   = ",                      input_data.Dissipation              )
    print( " Symmetry of system        = ",                      input_data.Symmetry                 )
    print( " The Runge-Kutta scheme in the time evolution   = ", input_data.Time_Evolution_Method    ) 
    print( " The finite-difference scheme in the simulation = ", input_data.Finite_Diffenence_Method )
    print(                                                                                           )
    print( " The static AMR grid type = ",                       input_data.static_grid_type         )
    print( " The moving AMR grid type = ",                       input_data.moving_grid_type         )
    print(                                                                                           )
    print( " The number of static AMR grid levels = ",           input_data.static_grid_level        )      
    print( " The number of moving AMR grid levels = ",           input_data.moving_grid_level        )      
    print( " The number of total  AMR grid levels = ",           input_data.grid_level               )
    print(                                                                                           )
    print( " The grid number of each static AMR grid level = ",  input_data.static_grid_number       )
    print( " The grid number of each moving AMR grid level = ",  input_data.moving_grid_number       )
    print(                                                                                           )
    print( " The scale for largest  static AMR grid in X direction = ", setup_data.maximal_domain_size_static_x )
    print( " The scale for largest  static AMR grid in Y direction = ", setup_data.maximal_domain_size_static_y )
    print( " The scale for largest  static AMR grid in Z direction = ", setup_data.maximal_domain_size_static_z )
    print( " The scale for smallest static AMR grid in X direction = ", setup_data.minimal_domain_size_static_x )
    print( " The scale for smallest static AMR grid in Y direction = ", setup_data.minimal_domain_size_static_y )
    print( " The scale for smallest static AMR grid in Z direction = ", setup_data.minimal_domain_size_static_z )
    print(                                                                                                      )

    
    if ( input_data.moving_grid_level > 0):
        print( " The scale for largest  moving AMR grid = ",     setup_data.maximal_domain_size_moving )
        print( " The scale for smallest moving AMR grid = ",     setup_data.minimal_domain_size_moving )

    print(                                                                                             )
    print( " The coarest resolution for static AMR grid = ",     setup_data.maximal_resolution_static  )
    print( " The finest  resolution for static AMR grid = ",     setup_data.minimal_resolution_static  )
    
    if ( input_data.moving_grid_level > 0):
        print( " The coarest resolution for moving AMR grid = ", setup_data.maximal_resolution_moving  )
        print( " The finest  resolution for moving AMR grid = ", setup_data.minimal_resolution_moving  )
    
    print(                                                                                                       )
    print( " The time refinement starts from AMR grid level = ", input_data.refinement_level                     )
    print( " The time interval in each step for coarest AMR grid during time evaluation = ", setup_data.TimeStep )
    print(                                                                                                       )

    if input_data.basic_grid_set == "Shell-Patch":
        print( " The Shell-Patch AMR grid structure is used in this simulation "                                 )
        print( " Shell-Patch grid number = ",                    input_data.shell_grid_number                    )
        print( " Shell-Patch grid minimal radius = ",            setup_data.minimal_domain_size_shellpatch       )
        print( " Shell-Patch grid maximal radius = ",            setup_data.maximal_domain_size_shellpatch       )
        print( " Shell-Patch grid radial  resolution = ",        setup_data.shellpatch_resolution_R              )
        print( " Shell-Patch grid angular resolution = ",        setup_data.shellpatch_resolution_phi,  \
                                                                 setup_data.shellpatch_resolution_theta          )
        print(                                                                                                   )
    elif input_data.basic_grid_set == "Patch":   
        print( " This simulation only uses the Patch AMR grid structure, the Shell-Patch is not used "  )
        print(                                                                                          )
    else:
        print( " The AMR grid structure setting is wrong !!! "                                          )
        print(                                                                                          )

    print( "------------------------------------------------------------------------------------------" ) 
    
    ## 文件输出
    
    filepath = os.path.join( File_directory, "AMSS_NCKU_resolution" )
    file0    = open(filepath, 'w')
         
    ## 下面输出英文信息
    
    print(                                                                                              file=file0 )
    print( " Printing the basic parameter and setting in the AMSS-NCKU simulation ",                    file=file0 )
    print(                                                                                              file=file0 )
    print( " The number of MPI processes in the AMSS-NCKU simulation = ", input_data.MPI_processes,     file=file0 ) 
    print(                                                                                              file=file0 )
    print( " The form of computational equation  = ",            input_data.Equation_Class,             file=file0 )
    print( " The initial data in this simulation = ",            input_data.Initial_Data_Method,        file=file0 )
    print(                                                                                              file=file0 )
    print( " Starting evolution time   = ",                      input_data.Start_Evolution_Time,       file=file0 ) 
    print( " Final evolution time      = ",                      input_data.Final_Evolution_Time,       file=file0 ) 
    print( " Maximal iteration number  = ",                      input_data.Evolution_Step_Number,      file=file0 )
    print( " Courant factor            = ",                      input_data.Courant_Factor,             file=file0 )
    print( " Strength of dissipation   = ",                      input_data.Dissipation,                file=file0 )
    print( " Symmetry of system        = ",                      input_data.Symmetry,                   file=file0 )
    print( " The Runge-Kutta scheme in the time evolution   = ", input_data.Time_Evolution_Method,      file=file0 ) 
    print( " The finite-difference scheme in the simulation = ", input_data.Finite_Diffenence_Method,   file=file0 )
    print(                                                                                              file=file0 )
    print( " The static AMR grid type = ",                       input_data.static_grid_type,           file=file0 )
    print( " The moving AMR grid type = ",                       input_data.moving_grid_type,           file=file0 )
    print(                                                                                              file=file0 )
    print( " The number of static AMR grid levels = ",           input_data.static_grid_level,          file=file0 )      
    print( " The number of moving AMR grid levels = ",           input_data.moving_grid_level,          file=file0 )      
    print( " The number of total  AMR grid levels = ",           input_data.grid_level,                 file=file0 )
    print(                                                                                              file=file0 )
    print( " The grid number of each static AMR grid level = ",  input_data.static_grid_number,         file=file0 )
    print( " The grid number of each moving AMR grid level = ",  input_data.moving_grid_number,         file=file0 )
    print(                                                                                              file=file0 )
    print( " The scale for largest  static AMR grid in X direction = ", setup_data.maximal_domain_size_static_x, file=file0 )
    print( " The scale for largest  static AMR grid in Y direction = ", setup_data.maximal_domain_size_static_y, file=file0 )
    print( " The scale for largest  static AMR grid in Z direction = ", setup_data.maximal_domain_size_static_z, file=file0 )
    print( " The scale for smallest static AMR grid in X direction = ", setup_data.minimal_domain_size_static_x, file=file0 )
    print( " The scale for smallest static AMR grid in Y direction = ", setup_data.minimal_domain_size_static_y, file=file0 )
    print( " The scale for smallest static AMR grid in Z direction = ", setup_data.minimal_domain_size_static_z, file=file0 )
    print(                                                                                                                  )
    
    if ( input_data.moving_grid_level > 0):
        print( " The scale for largest  moving AMR grid = ",     setup_data.maximal_domain_size_moving, file=file0 )
        print( " The scale for smallest moving AMR grid = ",     setup_data.minimal_domain_size_moving, file=file0 )

    print(                                                                                              file=file0 )
    print( " The coarest resolution for static AMR grid = ",     setup_data.maximal_resolution_static,  file=file0 )
    print( " The finest  resolution for static AMR grid = ",     setup_data.minimal_resolution_static,  file=file0 )
    
    if ( input_data.moving_grid_level > 0):
        print( " The coarest resolution for moving AMR grid = ", setup_data.maximal_resolution_moving,  file=file0 )
        print( " The finest  resolution for moving AMR grid = ", setup_data.minimal_resolution_moving,  file=file0 )
    
    print(                                                                                                        file=file0 )
    print( " The time refinement starts from AMR grid level = ", input_data.refinement_level,                     file=file0 )
    print( " The time interval in each step for coarest AMR grid during time evaluation = ", setup_data.TimeStep, file=file0 )
    print(                                                                                                        file=file0 )

    if input_data.basic_grid_set == "Shell-Patch":
        print( " The Shell-Patch AMR grid structure is used in this simulation ",                           file=file0 )
        print( " Shell-Patch grid number = ",                    input_data.shell_grid_number,              file=file0 )
        print( " Shell-Patch grid minimal radius = ",            setup_data.minimal_domain_size_shellpatch, file=file0 )
        print( " Shell-Patch grid maximal radius = ",            setup_data.maximal_domain_size_shellpatch, file=file0 )
        print( " Shell-Patch grid radial  resolution = ",        setup_data.shellpatch_resolution_R,        file=file0 )
        print( " Shell-Patch grid angular resolution = ",        setup_data.shellpatch_resolution_phi,  \
                                                                 setup_data.shellpatch_resolution_theta,    file=file0 )
        print(                                                                                              file=file0 )
    elif input_data.basic_grid_set == "Patch":   
        print( " This simulation only uses the Patch AMR grid structure, the Shell-Patch is not used ", file=file0 )
        print(                                                                                          file=file0 )
    else:
        print( " The AMR grid structure setting is wrong !!! ",                                         file=file0 )
        print(                                                                                          file=file0 )
        
    return

##################################################################




##################################################################

## 将以上基本信息写入 AMSS-NCKU 程序的输入文件

def generate_AMSSNCKU_input( input_language ): 

    ## 屏幕输出： 开始生成 AMSS-NCKU 输入文件
    if ( input_language == "Chinese" ):
        print(                                     )
        print( " 正在生成 AMSS-NCKU 程序的输入文件 " ) 
        print(                                     )
    elif ( input_language == "English" ): 
        print(                                                                                   )
        print( " Automatically generating the input parfile for AMSS-NCKU C++ program ABE.EXE. " ) 
        print(                                                                                   )

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU.input"), "w" ) 
    ## file1 = open( "AMSS-NCKU.input", "w" )  

    ## 输出 ABE 的相关设定

    print( file=file1 )
    print( "ABE::checkrun  =  0",                                          file=file1 )
    print( "ABE::checkfile =  bssn.chk",                                   file=file1 )
    print( "ABE::Steps     = ",         input_data.Evolution_Step_Number,  file=file1 )
    print( "ABE::StartTime = ",         input_data.Start_Evolution_Time,   file=file1 )
    print( "ABE::TotalTime = ",         input_data.Final_Evolution_Time,   file=file1 )
    print( "ABE::DumpTime  = ",         input_data.Dump_Time,              file=file1 )
    print( "ABE::d2DumpTime   = ",      input_data.D2_Dump_Time,           file=file1 )
    print( "ABE::CheckTime    = ",      input_data.Check_Time,             file=file1 )
    print( "ABE::AnalysisTime = ",      input_data.Analysis_Time,          file=file1 )
    print( "ABE::Courant      = ",      input_data.Courant_Factor,         file=file1 )

    if ( input_data.Symmetry == "octant-symmetry" ):
        print( "ABE::Symmetry     = 2 ",                                   file=file1 )
    elif ( input_data.Symmetry == "equatorial-symmetry" ):
        print( "ABE::Symmetry     = 1 ",                                   file=file1 )
    elif ( input_data.Symmetry == "no-symmetry" ):
        print( "ABE::Symmetry     = 0 ",                                   file=file1 )
    else :
        print( " Symmetry Setting Error" )

    print( "ABE::small dissipation = ",        input_data.Dissipation,     file=file1 )
    print( "ABE::big dissipation   = ",        input_data.Dissipation,     file=file1 )
    print( "ABE::shell dissipation = ",        input_data.Dissipation,     file=file1 )
    print( "ABE::Analysis Level    = ",        input_data.analysis_level,  file=file1 )
    print( "ABE::Max mode l        = ",        input_data.GW_L_max,        file=file1 )
    print( "ABE::detector number   = ",        input_data.Detector_Number, file=file1 )
    print( "ABE::farest detector position = ", input_data.Detector_Rmax,   file=file1 ) 
    print( f"ABE::detector distance = { (input_data.Detector_Rmax-input_data.Detector_Rmin) / (input_data.Detector_Number-1) }", \
           file=file1 )
    print( "ABE::cpu part   = ", input_data.CPU_Part,                      file=file1 ) 
    print( "ABE::gpu part   = ", input_data.GPU_Part,                      file=file1 )     
    print( "ABE::output dir = ", input_data.Output_directory,              file=file1 ) 
    
    if ( input_data.Initial_Data_Method == "Cao-Analytical" ):
        print( "ABE::ID Type    = -3",  file=file1 )
    elif ( input_data.Initial_Data_Method == "KerrSchild-Analytical" ):
        print( "ABE::ID Type    = -2", file=file1 )
    elif ( input_data.Initial_Data_Method == "Lousto-Analytical" ):
        print( "ABE::ID Type    = -1", file=file1 )
    elif ( input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
        print( "ABE::ID Type    = 0",  file=file1 )
    elif ( input_data.Initial_Data_Method == "Pablo-Olliptic" ):
        print( "ABE::ID Type    = 1",  file=file1 )
    else :
        print( " Initial Data Setting Error" )
        
    print( file=file1 )
    
    ## 输出 AHF 的相关设定
    
    print( "AHF::AHfindevery = ", input_data.AHF_Find_Every, file=file1 ) 
    print( "AHF::AHdumptime  = ", input_data.AHF_Dump_Time,  file=file1 ) 
    print(                                                   file=file1 )
    
    ## 输出其它设定
    print(                                                                                              file=file1 )
    print( "SurfaceIntegral::number of points for quarter sphere = ", input_data.quarter_sphere_number, file=file1 )
    print(                                                                                              file=file1 )
    
    ## 输出标量-张量-F(R)理论的设定
    ## 加入条件判断后，即使输入文件 AMSS_NCKU_Input.py 中没有  FR 相关的设定，也不会报错
    if (input_data.Equation_Class == "BSSN-EScalar"):
        ## 保留一定小数输出
        print( "FR::a2     = ", format(input_data.FR_a2,     '.2f'), file=file1 )
        print( "FR::l2     = ", format(input_data.FR_l2,     '.2f'), file=file1 )
        print( "FR::phi0   = ", format(input_data.FR_phi0,   '.8f'), file=file1 )
        print( "FR::r0     = ", format(input_data.FR_r0,     '.2f'), file=file1 )
        print( "FR::sigma0 = ", format(input_data.FR_sigma0, '.2f'), file=file1 )
        print(                                                       file=file1 )
    else:
        print(                                                       file=file1 )  
         
    file1.close()

    return 
    
##################################################################

        
