
#########################################################################
##
## 这个文件包含了数值相对论所需要的格点 | This file contains the grid points needed for numerical relativity
## 小曲 | Xiao Qu
## 2024/03/20
## 2026/02/27 修改 | Modified on 2026/02/27 
##
#########################################################################

import numpy                              ## 导入 numpy 包 | Import the numpy package 
import matplotlib.pyplot as plt           ## 导入 matplotlib 包中的 pyplot 函数，并用 plt 名称代表 matplotlib.pyplot | Import the pyplot module of the matplotlib package, with plt standing for matplotlib.pyplot
import os                                 ## 导入 os 包进行系统操作 | Import the os package for system operations

import AMSS_NCKU_Input   as input_data    ## 导入程序输入文件 | Import the program input file

import puncture_initialize                 ## 该文件包含 Puncture 的设定 | This file contains the Puncture setup

#########################################################################




#########################################################################

## 该函数用于处理数值计算的网格 | This function handles the grids for numerical computation

#########################################################################

def generate_initial_AMR_grid( input_language, puncture_data_input ):

    #################################################
    
    ## 导入 puncture 的信息 | Import the puncture information
    
    ## 为了更好的结构化，将这一步变为函数的输入 | For better structuring, this step has been made an input of the function
    ## puncture_data_input = puncture_initialize.generate_puncture_input_data( input_language )
    
    puncture = puncture_data_input.position_BH

    #################################################

    ## 输出网格的信息 | Output the grid information

    if ( input_language == "Chinese" ):
        print(                )   
        print( " 输出网格信息 " )
        print(                )
    elif ( input_language == "English" ):   
        print(                                     )   
        print( " Wirte Down The Grid Information " )
        print(                                     )
    
    if ( input_language == "Chinese" ):
        print( " 总网格层数   = ",   input_data.grid_level        )      ## 输出网格的总层数 | Output the total number of grid levels
        print( " 固定网格层数 = ",   input_data.static_grid_level )      ## 输出静态网格的层数 | Output the number of static grid levels
        print( " 可移动网格层数 = ", input_data.moving_grid_level )      ## 输出可移动网格的层数 | Output the number of moving grid levels
        print(                                                  )
    elif ( input_language == "English" ):     
        print( " Number of Total Grid Level = ",          input_data.grid_level        )      ## 输出网格的总层数 | Output the total number of grid levels
        print( " Number of Static Grid Level = ",         input_data.static_grid_level )      ## 输出静态网格的层数 | Output the number of static grid levels
        print( " Number of Moving Grid Level = ",         input_data.moving_grid_level )      ## 输出可移动网格的层数 | Output the number of moving grid levels
        ## print( " Number of Points in Each Grid Level = ", input_data.grid_number    )      ## 输出每一层每一维数的网格格点数目 | Output the number of grid points in each dimension of each grid level
        print(                                                                         )


    if ( input_language == "Chinese" ):
        print(                   )
        print( " 正在设定计算网格 " )
        print(                   )
    elif ( input_language == "English" ):
        print(                                         )
        print( " Setting the demanded numerical grid " )
        print(                                         )

    
    #################################################

    ## 初始化网格信息 | Initialize the grid information

    ## 初始化网格坐标最小值、最大值、格点数目 | Initialize the grid coordinate minima, maxima, and the number of grid points

    Grid_X_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 X 坐标的最小值 | Define a real array as the minimum X coordinate of each grid level
    Grid_X_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 X 坐标的最大值 | Define a real array as the maximum X coordinate of each grid level
    Grid_Y_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Y 坐标的最小值 | Define a real array as the minimum Y coordinate of each grid level
    Grid_Y_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Y 坐标的最大值 | Define a real array as the maximum Y coordinate of each grid level
    Grid_Z_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Z 坐标的最小值 | Define a real array as the minimum Z coordinate of each grid level
    Grid_Z_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Z 坐标的最大值 | Define a real array as the maximum Z coordinate of each grid level

    Grid_Resolution = numpy.zeros( input_data.grid_level ) # 定义实数数组，作为每层网格分辨率 | Define a real array as the resolution of each grid level

    largest_box_X_Max = input_data.largest_box_xyz_max[0] 
    largest_box_Y_Max = input_data.largest_box_xyz_max[1]
    largest_box_Z_Max = input_data.largest_box_xyz_max[2]
    largest_box_X_Min = input_data.largest_box_xyz_min[0] 
    largest_box_Y_Min = input_data.largest_box_xyz_min[1]
    largest_box_Z_Min = input_data.largest_box_xyz_min[2]

    ## 定义整数数组，作为每一层固定网格每一维的格点数 | Define an integer array as the number of grid points in each dimension of each static grid level
    static_grid_number_x = input_data.static_grid_number 
    static_grid_number_y = int( (largest_box_Y_Max - largest_box_Y_Min) * ( static_grid_number_x / (largest_box_X_Max-largest_box_X_Min) ) )
    static_grid_number_z = int( (largest_box_Z_Max - largest_box_Z_Min) * ( static_grid_number_x / (largest_box_X_Max-largest_box_X_Min) ) )
    
    ## 定义整数数组，作为每一层可移动网格每一维的格点数 | Define an integer array as the number of grid points in each dimension of each moving grid level
    moving_grid_number   = input_data.moving_grid_number   

    #################################################

    ## 初始化固定网格 | Initialize the static grids

    ## 将每个方向的网格数目自动调整为偶数 | Automatically adjust the number of grid points in each direction to be even
    ## Python 中使用 % 取余数 | In Python, % gives the remainder
    ## print(static_grid_number_x % 2)
    if ( (static_grid_number_x % 2) != 0) :
        static_grid_number_x = static_grid_number_x + 1
    if ( (static_grid_number_y % 2) != 0) :
        static_grid_number_y = static_grid_number_y + 1
    if ( (static_grid_number_z % 2) != 0) :
        static_grid_number_z = static_grid_number_z + 1
    ## 为了防止移动网格移动时边界无法与固定网格对其，我们进一步要求每个方向的网格数目自动调整为 4 的倍数 | To prevent the moving-grid boundaries from failing to align with the static grids when the grids move, we further require the number of grid points in each direction to be adjusted to a multiple of 4
    if ( (static_grid_number_x % 4) != 0) :
        static_grid_number_x = static_grid_number_x + 2
    if ( (static_grid_number_y % 4) != 0) :
        static_grid_number_y = static_grid_number_y + 2
    if ( (static_grid_number_z % 4) != 0) :
        static_grid_number_z = static_grid_number_z + 2
    '''
    ## 为了防止移动网格移动时边界无法与固定网格对其，我们进一步要求每个方向的网格数目自动调整为 8 的倍数 | To prevent the moving-grid boundaries from failing to align with the static grids when the grids move, we further require the number of grid points in each direction to be adjusted to a multiple of 8
    if ( (static_grid_number_x % 8) != 0) :
        static_grid_number_x = static_grid_number_x + 4
    if ( (static_grid_number_y % 8) != 0) :
        static_grid_number_y = static_grid_number_y + 4
    if ( (static_grid_number_z % 8) != 0) :
        static_grid_number_z = static_grid_number_z + 4
    '''

    ##  定义实数数组，维度 grid_number * static_grid_level，分别作为每一层固定网格的 X Y Z 坐标 | Define real arrays of dimension grid_number * static_grid_level as the X, Y, Z coordinates of each static grid level
    Static_Grid_X = numpy.zeros( (input_data.static_grid_level, static_grid_number_x+1) )   
    Static_Grid_Y = numpy.zeros( (input_data.static_grid_level, static_grid_number_y+1) )   
    Static_Grid_Z = numpy.zeros( (input_data.static_grid_level, static_grid_number_z+1) )  

    #################################################

    ## 初始化可移动网格 | Initialize the moving grids

    ##  定义实数数组，维度 grid_number * puncture_number * moving_grid_level，作为每一层可移动网格的 X Y Z 坐标 | Define real arrays of dimension grid_number * puncture_number * moving_grid_level as the X, Y, Z coordinates of each moving grid level
    Moving_Grid_X = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) )
    Moving_Grid_Y = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) ) 
    Moving_Grid_Z = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) )

    #################################################

    ## 初始化可移动网格的坐标最小值、最大值 | Initialize the coordinate minima and maxima of the moving grids

    Moving_Grid_X_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) ) 
    Moving_Grid_X_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
    Moving_Grid_Y_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
    Moving_Grid_Y_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
    Moving_Grid_Z_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )
    Moving_Grid_Z_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )          

    #################################################

    ## 设置每层网格的网格分辨率 | Set the resolution of each grid level

    for i in range(input_data.static_grid_level) :
        if i==0:
            Grid_Resolution[i] = ( largest_box_X_Max - largest_box_X_Min ) / static_grid_number_x
        else:
            Grid_Resolution[i] = Grid_Resolution[i-1] / input_data.devide_factor
    
    for j in range(input_data.moving_grid_level) : 
        i = j + input_data.static_grid_level
        Grid_Resolution[i] = Grid_Resolution[i-1] / input_data.devide_factor

    #################################################

    ## 根据输入文件，设置每层固定网格坐标的最小值和最大值 | According to the input file, set the coordinate minima and maxima of each static grid level

    ## 设置第一层固定网格的坐标最小值与最大值 | Set the coordinate minima and maxima of the first static grid level
    Grid_X_Min[0] = largest_box_X_Min 
    Grid_X_Max[0] = largest_box_X_Max
    Grid_Y_Min[0] = largest_box_Y_Min  
    Grid_Y_Max[0] = largest_box_Y_Max
    Grid_Z_Min[0] = largest_box_Z_Min
    Grid_Z_Max[0] = largest_box_Z_Max
    ## 重新调整 yx 的最大坐标，确保 xyz 三个方向的分辨率一致 | Readjust the maximum yz coordinates to ensure the same resolution in the x, y and z directions
    Grid_Y_Max[0] = Grid_Y_Min[0] + Grid_Resolution[0] * static_grid_number_y
    Grid_Z_Max[0] = Grid_Z_Min[0] + Grid_Resolution[0] * static_grid_number_z
    ## 如果系统存在对称性，还要进行进一步的调整 | If the system has symmetries, further adjustments are still needed
    if ( input_data.Symmetry == "equatorial-symmetry" ):
        Grid_Z_Min[0] = - Grid_Resolution[0] * static_grid_number_z / 2
        Grid_Z_Max[0] = + Grid_Resolution[0] * static_grid_number_z / 2
    elif ( input_data.Symmetry == "octant-symmetry" ):
        Grid_X_Min[0] = - Grid_Resolution[0] * static_grid_number_x / 2
        Grid_X_Max[0] = + Grid_Resolution[0] * static_grid_number_x / 2
        Grid_Y_Min[0] = - Grid_Resolution[0] * static_grid_number_y / 2
        Grid_Y_Max[0] = + Grid_Resolution[0] * static_grid_number_y / 2
        Grid_Z_Min[0] = - Grid_Resolution[0] * static_grid_number_z / 2
        Grid_Z_Max[0] = + Grid_Resolution[0] * static_grid_number_z / 2

    ## print( " Grid_Y_Max[0] = ", Grid_Y_Max[0] )

    if ( input_language == "Chinese" ):
        print(                                                   )
        print( " 正在调整固定网格格点的位置，使得坐标原点位于固定网格上 " )
        print(                                                   )
    elif ( input_language == "English" ):
        print(                                                                                                   )
        print( " adjusting the static gird points, making the original point (0,0,0) to the static gird points " )
        print(                                                                                                   )

    ## 设置其它层固定网格的坐标最大值与最小值 | Set the coordinate maxima and minima of the other static grid levels
    for i in range(input_data.static_grid_level-1) :
        
        ## 如果坐标原点不在最外层固定网格上，调整最外层固定网格，使得原点位于最外层网格上 | If the coordinate origin is not on the outermost static grid, adjust the outermost static grid so that the origin lies on it
        if i==0:
            for nn in range(static_grid_number_x):
                if (Grid_X_Min[i] + nn*Grid_Resolution[i]) < 0.0 < (Grid_X_Min[i] + (nn+1)*Grid_Resolution[i]):
                    print( " before adjust: Grid X_min = ", Grid_X_Min[i] )
                    print( " before adjust: Grid X_max = ", Grid_X_Max[i] )
                    grid_adjust   = Grid_X_Min[i] + (nn+1)*Grid_Resolution[i]
                    Grid_X_Min[i] = Grid_X_Min[i] - grid_adjust
                    Grid_X_Max[i] = Grid_X_Max[i] - grid_adjust
                    print( " after adjust: Grid X_min = ", Grid_X_Min[i] )
                    print( " after adjust: Grid X_max = ", Grid_X_Max[i] )
            for nn in range(static_grid_number_y):
                if (Grid_Y_Min[i] + nn*Grid_Resolution[i]) < 0.0 < (Grid_Y_Min[i] + (nn+1)*Grid_Resolution[i]):
                    print( " before adjust: Grid Y_min = ", Grid_Y_Min[i] )
                    print( " before adjust: Grid Y_max = ", Grid_Y_Max[i] )
                    grid_adjust   = Grid_Y_Min[i] + (nn+1)*Grid_Resolution[i]
                    Grid_Y_Min[i] = Grid_Y_Min[i] - grid_adjust
                    Grid_Y_Max[i] = Grid_Y_Max[i] - grid_adjust
                    print( " after adjust: Grid Y_min = ", Grid_Y_Min[i] )
                    print( " after adjust: Grid Y_max = ", Grid_Y_Max[i] )
            for nn in range(static_grid_number_z):
                if (Grid_Z_Min[i] + nn*Grid_Resolution[i]) < 0.0 < (Grid_Z_Min[i] + (nn+1)*Grid_Resolution[i]):
                    print( " before adjust: Grid Z_min = ", Grid_Z_Min[i] )
                    print( " before adjust: Grid Z_max = ", Grid_Z_Max[i] )
                    grid_adjust   = Grid_X_Min[i] + (nn+1)*Grid_Resolution[i]
                    Grid_Z_Min[i] = Grid_Z_Min[i] - grid_adjust
                    Grid_Z_Max[i] = Grid_Z_Max[i] - grid_adjust
                    print( " after adjust: Grid Z_min = ", Grid_Z_Min[i] )
                    print( " after adjust: Grid Z_max = ", Grid_Z_Max[i] )
        
        ## 每层固定网格最小值（或最大值）为上一层固定网格最小值（或最大值）除以因子 devide_factor | The minimum (or maximum) of each static grid level equals that of the previous level divided by the factor devide_factor
        Grid_X_Min[i+1] = Grid_X_Min[i] / input_data.devide_factor    
        Grid_X_Max[i+1] = Grid_X_Max[i] / input_data.devide_factor    
        Grid_Y_Min[i+1] = Grid_Y_Min[i] / input_data.devide_factor
        Grid_Y_Max[i+1] = Grid_Y_Max[i] / input_data.devide_factor
        Grid_Z_Min[i+1] = Grid_Z_Min[i] / input_data.devide_factor
        Grid_Z_Max[i+1] = Grid_Z_Max[i] / input_data.devide_factor
    

    ## 固定网格与可移动网格的格点数目不同，为保证相同分辨率，需要引入网格大小调节因子 | The static and moving grids have different numbers of grid points; to guarantee the same resolution, a grid-size adjustment factor is introduced
    adjust_factor = input_data.moving_grid_number / input_data.static_grid_number

    ## 设置第一层可移动网格的坐标最大值与最小值 | Set the coordinate maxima and minima of the first moving grid level
    i = input_data.static_grid_level 
    if (i < input_data.grid_level):
 
        ## 原来的设定，移动网格可能变为长方体网格 | The old setup, in which the moving grid could become a cuboid grid
        ## Grid_X_Min[i] = ( Grid_X_Min[i-1] / input_data.devide_factor ) * adjust_factor 
        ## Grid_X_Max[i] = ( Grid_X_Max[i-1] / input_data.devide_factor ) * adjust_factor 
        ## Grid_Y_Min[i] = ( Grid_Y_Min[i-1] / input_data.devide_factor ) * adjust_factor
        ## Grid_Y_Max[i] = ( Grid_Y_Max[i-1] / input_data.devide_factor ) * adjust_factor
        ## Grid_Z_Min[i] = ( Grid_Z_Min[i-1] / input_data.devide_factor ) * adjust_factor
        ## Grid_Z_Max[i] = ( Grid_Z_Max[i-1] / input_data.devide_factor ) * adjust_factor

        ## 现在的设定，总是保证移动网格为立方网格 | The current setup, which always keeps the moving grid cubic
        Grid_X_Min[i] = - ( 0.5 * (Grid_X_Max[i-1] - Grid_X_Min[i-1]) / input_data.devide_factor ) * adjust_factor
        Grid_X_Max[i] = - Grid_X_Min[i] 

        Grid_Y_Min[i] = Grid_X_Min[i]
        Grid_Y_Max[i] = Grid_X_Max[i]
        Grid_Z_Min[i] = Grid_X_Min[i]
        Grid_Z_Max[i] = Grid_X_Max[i]

        ## print( " Grid_X_Max[i] = ", Grid_X_Max[i] )
        ## print( " Grid_Y_Max[i] = ", Grid_Y_Max[i] )

    ## 设置其它层可移动网格的坐标最大值与最小值 | Set the coordinate maxima and minima of the other moving grid levels
    for j in range(input_data.moving_grid_level-1) :
        k = input_data.static_grid_level + j
        Grid_X_Min[k+1] = Grid_X_Min[k] / input_data.devide_factor    
        Grid_X_Max[k+1] = Grid_X_Max[k] / input_data.devide_factor    
        Grid_Y_Min[k+1] = Grid_Y_Min[k] / input_data.devide_factor
        Grid_Y_Max[k+1] = Grid_Y_Max[k] / input_data.devide_factor
        Grid_Z_Min[k+1] = Grid_Z_Min[k] / input_data.devide_factor
        Grid_Z_Max[k+1] = Grid_Z_Max[k] / input_data.devide_factor
    
    ## 设置最外层球壳格点 shell patch 的坐标最大值和最小值 | Set the coordinate maxima and minima of the outermost shell-patch grid points

    Shell_R_Resolution = Grid_Resolution[0]
    Shell_R_Min        = largest_box_X_Max
    Shell_R_Max        = Shell_R_Min + Grid_Resolution[0] * input_data.shell_grid_number[2]

    #################################################


    #################################################

    ## 设置每层网格的坐标值 | Set the coordinate values of each grid level

    #################################################

    ## 设置固定网格的格点坐标坐标 | Set the grid-point coordinates of the static grids

    ## 线性格点 | Linear grid points

    if input_data.static_grid_type == 'Linear' : 

        for i in range(input_data.static_grid_level):
            Static_Grid_X[i] = numpy.linspace( Grid_X_Min[i], Grid_X_Max[i], static_grid_number_x+1 ) 
            Static_Grid_Y[i] = numpy.linspace( Grid_Y_Min[i], Grid_Y_Max[i], static_grid_number_y+1 ) 
            Static_Grid_Z[i] = numpy.linspace( Grid_Z_Min[i], Grid_Z_Max[i], static_grid_number_z+1 )
         ## 利用 numpy 设置线性格点，参数表示起始值 Rmin，最大值 Rmax，格点数目 Rnum | Use numpy to set up linear grid points, where the parameters denote the starting value Rmin, the maximum value Rmax, and the number of grid points Rnum
         ## 注意，如果是线性格点，则格点坐标最大值为 GridMax；如果是对数格点，则格点坐标最大值为 e^{GridMax} | Note: for linear grid points the maximum grid coordinate is GridMax, while for logarithmic grid points it is e^{GridMax}

    else:
        if ( input_language == "Chinese" ):
            print(                                )
            print( " 固定网格：网格类型设置错误 ！！ " )
            print(                                )
        elif ( input_language == "English" ):
            print(                                                       )
            print( " Static Grid Error: Grid Type is Undifined !!!!!!! " )
            print(                                                       )

    #################################################
  
    ## 设置可移动网格的格点坐标坐标 | Set the grid-point coordinates of the moving grids  

    if ( input_language == "Chinese" ):
        print(                                                                     )
        print( " 正在调整移动网格中心的位置，以便同时与 Puncture 位置以及上一层粗网格相匹配 " )
        print(                                                                     )
    elif ( input_language == "English" ):
        print(                                                                                                            )
        print( " adjusting the moving gird points, ensuring the alliance of moving grids points and static grids points " )
        print(                                                                                                            )

    ## 为了让每一层可移动网格的边界与上一层网格重合，要对 Puncture 的位置进行调整 | To make the boundary of each moving grid level coincide with the previous grid level, the Puncture position is adjusted
    adjust_puncture = numpy.zeros( (input_data.puncture_number, 3) )

    ## 线性格点 | Linear grid points
  
    if ( input_data.moving_grid_type == "Linear" ): 
    
        ## 对可移动网格层数进行循环 | Loop over the moving grid levels
        for j in range(input_data.moving_grid_level) : 

            i = j + input_data.static_grid_level
        
            ## 对 Puncture 数目进行循环 | Loop over the number of Punctures
            for k in range(input_data.puncture_number) :
            
                ## 将 Puncture 的位置进行调整 | Adjust the position of the Puncture
                if j==0 :
            
                    level0 = input_data.static_grid_level - 1  ## 引入新变量，避免代码过长 | Introduce a new variable to avoid overly long lines
                
                    for m in range(static_grid_number_x) : 
                        if ( Static_Grid_X[level0, m] <= puncture[k,0] <= Static_Grid_X[level0, m+1] ):
                            if ( abs( puncture[k,0] - Static_Grid_X[level0, m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,0] = Static_Grid_X[level0, m]
                            elif ( abs( puncture[k,0] - Static_Grid_X[level0, m+1 ] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,0] = Static_Grid_X[level0, m+1 ]
                            else:
                                adjust_puncture[k,0] = ( Static_Grid_X[level0, m] + Static_Grid_X[level0, m+1] ) / 2.0
                
                    for m in range(static_grid_number_y) :
                        if ( Static_Grid_Y[level0, m] <= puncture[k,1] <= Static_Grid_Y[level0, m+1] ):
                            if ( abs( puncture[k,1] - Static_Grid_Y[level0, m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,1] = Static_Grid_Y[level0, m]
                            elif ( abs( puncture[k,1] - Static_Grid_Y[level0, m+1] ) < ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,1] = Static_Grid_Y[level0, m+1]
                            else:
                                adjust_puncture[k,1] = ( Static_Grid_Y[level0, m] + Static_Grid_Y[level0, m+1] ) / 2.0
                
                    for m in range(static_grid_number_z) :
                        if ( Static_Grid_Z[level0, m] <= puncture[k,2] <= Static_Grid_Z[level0, m+1] ):
                            if ( abs( puncture[k,2] - Static_Grid_Z[level0, m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,2] = Static_Grid_Z[level0, m]
                            elif ( abs( puncture[k,2] - Static_Grid_Z[level0, m+1] ) < ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,2] = Static_Grid_Z[level0, m+1]
                            else:
                                adjust_puncture[k,2] = ( Static_Grid_Z[level0, m] + Static_Grid_Z[level0, m+1] ) / 2.0


                elif j>0 :
                    for m in range(moving_grid_number) :
                
                        if ( Moving_Grid_X[j-1,k,m] <= puncture[k,0] <= Moving_Grid_X[j-1,k,m+1] ):
                            if ( abs( puncture[k,0] - Moving_Grid_X[j-1,k,m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,0] = Moving_Grid_X[j-1,k,m]
                            elif ( abs( puncture[k,0] - Moving_Grid_X[j-1,k,m+1] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,0] = Moving_Grid_X[j-1,k,m+1]
                            else:
                                adjust_puncture[k,0] = ( Moving_Grid_X[j-1,k,m] + Moving_Grid_X[j-1,k,m+1] ) / 2.0
                
                        if ( Moving_Grid_Y[j-1,k,m] <= puncture[k,1] <= Moving_Grid_Y[j-1,k,m+1] ):
                            if ( abs( puncture[k,1] - Moving_Grid_Y[j-1,k,m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,1] = Moving_Grid_Y[j-1,k,m]
                            elif ( abs( puncture[k,1] - Moving_Grid_Y[j-1,k,m+1] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,1] = Moving_Grid_Y[j-1,k,m+1]
                            else:
                                adjust_puncture[k,1] = ( Moving_Grid_Y[j-1,k,m] + Moving_Grid_Y[j-1,k,m+1] ) / 2.0

                        if ( Moving_Grid_Z[j-1,k,m] <= puncture[k,2] <= Moving_Grid_Z[j-1,k,m+1] ):
                            if ( abs( puncture[k,2] - Moving_Grid_Z[j-1,k,m] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,2] = Moving_Grid_Z[j-1,k,m]
                            elif ( abs( puncture[k,2] - Moving_Grid_Z[j-1,k,m+1] )  <  ( Grid_Resolution[i]/2.0 ) ):
                                adjust_puncture[k,2] = Moving_Grid_Z[j-1,k,m+1]
                            else:
                                adjust_puncture[k,2] = ( Moving_Grid_Z[j-1,k,m] + Moving_Grid_Z[j-1,k,m+1] ) / 2.0

                else:
                    print( " Adjusting puncture position to compatable with coaser grid !  Error !!! " )
                ## Puncture 的位置调整完毕 | The adjustment of the Puncture position is finished
            
                ## 为了防止 C++ 程度读入数据出错 | To prevent the C++ program from reading in data incorrectly
                ## 如果遇到很小的数 1e-10，保留 2 位小数，不使用科学技术法 | If a very small number such as 1e-10 occurs, keep 2 decimal places and do not use scientific notation
                if ( abs(adjust_puncture[k,0]) < 1e-10 ):
                    adjust_puncture[k,0] = 0.00
                    ## adjust_puncture[k,0] = f"{ adjust_puncture[k,0]:.2f }
                if ( abs(adjust_puncture[k,1]) < 1e-10 ):
                    adjust_puncture[k,1] = 0.00
                    ## adjust_puncture[k,1] = f"{ adjust_puncture[k,1]:.2f }
                if ( abs(adjust_puncture[k,2]) < 1e-10 ):
                    adjust_puncture[k,2] = 0.00
                    ## adjust_puncture[k,2] = f"{ adjust_puncture[k,2]:.2f }
            
                ## 第 j 层可移动网格的 XYZ 坐标最小值（或最大值） = 第 i 层网格的 XYZ 坐标最小值（或最大值） + 第 k 个 puncture 的 XYZ 坐标 | The XYZ coordinate minimum (or maximum) of the j-th moving grid level = the XYZ coordinate minimum (or maximum) of the i-th grid level + the XYZ coordinate of the k-th puncture
                Moving_Grid_X_Min[j,k] = adjust_puncture[k,0] + Grid_X_Min[i]  
                Moving_Grid_X_Max[j,k] = adjust_puncture[k,0] + Grid_X_Max[i]  
                Moving_Grid_Y_Min[j,k] = adjust_puncture[k,1] + Grid_Y_Min[i]  
                Moving_Grid_Y_Max[j,k] = adjust_puncture[k,1] + Grid_Y_Max[i]  
                Moving_Grid_Z_Min[j,k] = adjust_puncture[k,2] + Grid_Z_Min[i]  
                Moving_Grid_Z_Max[j,k] = adjust_puncture[k,2] + Grid_Z_Max[i]
            
                ## 为了防止 C++ 程度读入数据出错 | To prevent the C++ program from reading in data incorrectly
                ## 如果遇到很小的数 1e-10，保留 2 位小数，不使用科学技术法 | If a very small number such as 1e-10 occurs, keep 2 decimal places and do not use scientific notation
                if ( abs(Moving_Grid_X_Min[j,k]) < 1e-10 ):
                    Moving_Grid_X_Min[j,k] = 0.00 
                if ( abs(Moving_Grid_X_Max[j,k]) < 1e-10 ):
                    Moving_Grid_X_Max[j,k] = 0.00 
                
                if ( abs(Moving_Grid_Y_Min[j,k]) < 1e-10 ):
                    Moving_Grid_Y_Min[j,k] = 0.00 
                if ( abs(Moving_Grid_Y_Max[j,k]) < 1e-10 ):
                    Moving_Grid_Y_Max[j,k] = 0.00 
                
                if ( abs(Moving_Grid_Z_Min[j,k]) < 1e-10 ):
                    Moving_Grid_Z_Min[j,k] = 0.00 
                if ( abs(Moving_Grid_Z_Max[j,k]) < 1e-10 ):
                    Moving_Grid_Z_Max[j,k] = 0.00 
            
                if ( input_language == "Chinese" ):
                    print( f" 移动网格位置调整  adjust_puncture[{i},{k},0] = { adjust_puncture[k,0] } " )
                    print( f" 移动网格位置调整  adjust_puncture[{i},{k},1] = { adjust_puncture[k,1] } " )
                    print( f" 移动网格位置调整  adjust_puncture[{i},{k},2] = { adjust_puncture[k,2] } " )
                elif ( input_language == "English" ):
                    print( f" adjust_puncture[{i},{k},0] = { adjust_puncture[k,0] } " )
                    print( f" adjust_puncture[{i},{k},1] = { adjust_puncture[k,1] } " )
                    print( f" adjust_puncture[{i},{k},2] = { adjust_puncture[k,2] } " )
            
                ## 利用 numpy 设置线性格点，参数表示起始值 Rmin，最大值 Rmax，格点数目 Rnum | Use numpy to set up linear grid points, where the parameters denote the starting value Rmin, the maximum value Rmax, and the number of grid points Rnum
                Moving_Grid_X[j,k] = numpy.linspace( Moving_Grid_X_Min[j,k], Moving_Grid_X_Max[j,k], moving_grid_number + 1 )  
                Moving_Grid_Y[j,k] = numpy.linspace( Moving_Grid_Y_Min[j,k], Moving_Grid_Y_Max[j,k], moving_grid_number + 1 )
                Moving_Grid_Z[j,k] = numpy.linspace( Moving_Grid_Z_Min[j,k], Moving_Grid_Z_Max[j,k], moving_grid_number + 1 )
        
    else:
        if ( input_language == "Chinese" ):
            print(                                  )
            print( " 可移动网格：网格类型设置错误 ！！ " )
            print(                                  )
        elif ( input_language == "English" ):
            print(                                                       )
            print( " Moving Grid Error: Grid Type is Undifined !!!!!!! " )
            print(                                                       )

    if ( input_language == "Chinese" ):
        print(                          )
        print( " 移动网格中心位置调整完毕 " )
        print(                          )
    elif ( input_language == "English" ):
        print(                                                         )
        print( " adjustment of moving grid points has been completed " )
        print(                                                         )

    ## 测试 | Test
    ## print( Moving_Grid_X_Min[0,0] ) 

    #################################################
    
    ## 设定 class 作为结构体，用来返回变量 | Define a class as a struct to return the variables
    
    class Initial_AMR_Grid_Class:
        def __init__( self, Static_Grid_X, Static_Grid_Y, Static_Grid_Z, Moving_Grid_X,     Moving_Grid_Y,     Moving_Grid_Z,     \
                            Grid_X_Max,    Grid_Y_Max,    Grid_Z_Max,    Moving_Grid_X_Max, Moving_Grid_Y_Max, Moving_Grid_Z_Max, \
                            Grid_X_Min,    Grid_Y_Min,    Grid_Z_Min,    Moving_Grid_X_Min, Moving_Grid_Y_Min, Moving_Grid_Z_Min, \
                            static_grid_number_x,         static_grid_number_y,             static_grid_number_z,                 \
                            Shell_R_Min,   Shell_R_Max ):
        
            self.Shell_R_Min   = Shell_R_Min
            self.Shell_R_Max   = Shell_R_Max
            self.Static_Grid_X = Static_Grid_X
            self.Static_Grid_Y = Static_Grid_Y
            self.Static_Grid_Z = Static_Grid_Z
            self.Moving_Grid_X = Moving_Grid_X
            self.Moving_Grid_Y = Moving_Grid_Y
            self.Moving_Grid_Z = Moving_Grid_Z
            self.Grid_X_Max    = Grid_X_Max
            self.Grid_Y_Max    = Grid_Y_Max
            self.Grid_Z_Max    = Grid_Z_Max
            self.Grid_X_Min    = Grid_X_Min
            self.Grid_Y_Min    = Grid_Y_Min
            self.Grid_Z_Min    = Grid_Z_Min
            self.Moving_Grid_X_Max = Moving_Grid_X_Max
            self.Moving_Grid_Y_Max = Moving_Grid_Y_Max
            self.Moving_Grid_Z_Max = Moving_Grid_Z_Max
            self.Moving_Grid_X_Min = Moving_Grid_X_Min
            self.Moving_Grid_Y_Min = Moving_Grid_Y_Min
            self.Moving_Grid_Z_Min = Moving_Grid_Z_Min
            self.static_grid_number_x = static_grid_number_x
            self.static_grid_number_y = static_grid_number_y
            self.static_grid_number_z = static_grid_number_z

    
    ## 函数的新版返回不再需要这个 class | The new-style return of the function no longer needs this class
    '''
    class Puncture_Output_Class:
        def __init__( self, BBH_M1,                BBH_M2,          \
                            distance_d0,           ellipticity_e0,  \
                            dimensionless_mass_BH, charge_Q_BH,     \
                            position_BH,           momentum_BH,     \
                            angular_momentum_BH ):
            self.BBH_M1                = BBH_M1
            self.BBH_M2                = BBH_M2
            self.distance_d0           = distance_d0
            self.ellipticity_e0        = ellipticity_e0
            self.dimensionless_mass_BH = dimensionless_mass_BH
            self.charge_Q_BH           = charge_Q_BH
            self.position_BH           = position_BH
            self.momentum_BH           = momentum_BH
            self.angular_momentum_BH   = angular_momentum_BH
    '''
            
    #################################################
    
    '''
    BBH_M1                = puncture_data_input.BBH_M1
    BBH_M2                = puncture_data_input.BBH_M2
    distance_d0           = puncture_data_input.distance_d0
    ellipticity_e0        = puncture_data_input.ellipticity_e0
    dimensionless_mass_BH = puncture_data_input.dimensionless_mass_BH
    charge_Q_BH           = puncture_data_input.charge_Q_BH
    position_BH           = puncture_data_input.position_BH
    momentum_BH           = puncture_data_input.momentum_BH
    angular_momentum_BH   = puncture_data_input.angular_momentum_BH
            
    ## 创建 class 实例 | Create a class instance
    puncture_data_output = Puncture_Output_Class( BBH_M1,                BBH_M2,         \
                                                  distance_d0,           ellipticity_e0, \
                                                  dimensionless_mass_BH, charge_Q_BH,    \
                                                  position_BH,           momentum_BH,    \
                                                  angular_momentum_BH )
    '''
    
    ## 以下为测试，检查是否创建正确的 class 实例 | The following is a test to check whether the correct class instance has been created
    '''
    print( puncture_data_output.dimensionless_mass_BH      )
    print( puncture_data_output.charge_Q_BH )
    print( puncture_data_output.position_BH )
    print( puncture_data_output.momentum_BH )
    print( puncture_data_output.angular_momentum_BH )
    '''
            
    ## 创建 class 实例 | Create a class instance
    initial_AMR_grid_data = Initial_AMR_Grid_Class( Static_Grid_X, Static_Grid_Y, Static_Grid_Z, Moving_Grid_X,     Moving_Grid_Y,     Moving_Grid_Z,     \
                                                    Grid_X_Max,    Grid_Y_Max,    Grid_Z_Max,    Moving_Grid_X_Max, Moving_Grid_Y_Max, Moving_Grid_Z_Max, \
                                                    Grid_X_Min,    Grid_Y_Min,    Grid_Z_Min,    Moving_Grid_X_Min, Moving_Grid_Y_Min, Moving_Grid_Z_Min, \
                                                    static_grid_number_x,         static_grid_number_y,             static_grid_number_z,                 \
                                                    Shell_R_Min,   Shell_R_Max )
    
    ## 以下为测试，检查是否创建正确的 class 实例 | The following is a test to check whether the correct class instance has been created
    '''
    print( initial_AMR_grid_data.static_grid_number_x )
    print( initial_AMR_grid_data.static_grid_number_y )
    print( initial_AMR_grid_data.static_grid_number_z )
    print( initial_AMR_grid_data.Grid_X_Min )
    print( initial_AMR_grid_data.Grid_X_Max )
    print( initial_AMR_grid_data.Grid_Y_Min )
    print( initial_AMR_grid_data.Grid_Y_Max )
    print( initial_AMR_grid_data.Grid_Z_Min )
    print( initial_AMR_grid_data.Grid_Z_Max )
    print( initial_AMR_grid_data.Shell_R_Min )   
    print( initial_AMR_grid_data.Shell_R_Max )
    '''
    
    #######################################################
    
    ## 函数结束 | End of the function
    
    ## 旧版返回 | Old-style return
    ## return initial_AMR_grid_data, puncture_data_output
    
    ## 新版返回不需要 puncture_data | The new-style return does not need puncture_data
    return initial_AMR_grid_data

#########################################################################



#########################################################################

## 该函数画出初始的网格 | This function plots the initial grids
    
def plot_initial_grid( input_language, initial_AMR_grid_data ):

    ## 屏幕输出：开始画出初始网格 | Screen output: start plotting the initial grid
    if ( input_language == "Chinese" ):
        print(                        )
        print( " 正在画出初始格点的图像 " ) 
        print(                        )
    elif ( input_language == "English" ):
        print(                                                      )
        print( " Schematically plot the numerical grid structure. " ) 
        print(                                                      )

    ## 依次画出最后的网格 | Plot the final grids one by one

    if (input_data.static_grid_level > 0):
        X0, Y0 = numpy.meshgrid( initial_AMR_grid_data.Static_Grid_X[0], initial_AMR_grid_data.Static_Grid_Y[0] )
        plt.plot( X0, Y0,                         # 画出第 0 层静态格点 | Plot the level-0 static grid points
                  color='brown',  	          # 全部点设置为红棕色 | Set all points to reddish brown
                  marker='.',  	                  # 点的形状为圆点 | The marker shape is a round dot
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
              
    if (input_data.static_grid_level > 1):
        X1, Y1 = numpy.meshgrid( initial_AMR_grid_data.Static_Grid_X[1], initial_AMR_grid_data.Static_Grid_Y[1] )
        plt.plot( X1, Y1,                         # 画出第 1 层静态格点 | Plot the level-1 static grid points
                  color='red',                    # 全部点设置为红色 | Set all points to red
                  marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
                  
    if (input_data.static_grid_level > 2):
        X2, Y2 = numpy.meshgrid( initial_AMR_grid_data.Static_Grid_X[2], initial_AMR_grid_data.Static_Grid_Y[2] )
        plt.plot( X2, Y2,                         # 画出第 2 层静态格点 | Plot the level-2 static grid points
                  color='orange',                 # 全部点设置为橙色 | Set all points to orange
                  marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
                  
    if (input_data.static_grid_level > 3):
        X3, Y3 = numpy.meshgrid( initial_AMR_grid_data.Static_Grid_X[3], initial_AMR_grid_data.Static_Grid_Y[3] )
        plt.plot( X3, Y3,                         # 画出第 3 层静态格点 | Plot the level-3 static grid points
                  color='yellow',                 # 全部点设置为黄色 | Set all points to yellow
                  marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
                  
    if (input_data.static_grid_level > 4):
        X4, Y4 = numpy.meshgrid( initial_AMR_grid_data.Static_Grid_X[4], initial_AMR_grid_data.Static_Grid_Y[4] )
        plt.plot( X4, Y4,                         # 画出第 4 层静态格点 | Plot the level-4 static grid points
                  color='greenyellow',            # 全部点设置为黄绿色 | Set all points to yellow-green
                  marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
                  
    ## 依次画出可移动网格 | Plot the moving grids one by one

    if (input_data.moving_grid_level > 0):
        for k in range(input_data.puncture_number):
            Xk0, Yk0 = numpy.meshgrid( initial_AMR_grid_data.Moving_Grid_X[0,k], initial_AMR_grid_data.Moving_Grid_Y[0,k] )
            plt.plot( Xk0, Yk0,                       # 画出第 0 层可移动格点 | Plot the level-0 moving grid points
                      color='cyan',  	              # 全部点设置为青蓝色 | Set all points to cyan
                      marker='.',  	              # 点的形状为圆点 | The marker shape is a round dot
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
 
    if (input_data.moving_grid_level > 1):
        for k in range(input_data.puncture_number):
            Xk1, Yk1 = numpy.meshgrid( initial_AMR_grid_data.Moving_Grid_X[1,k], initial_AMR_grid_data.Moving_Grid_Y[1,k] )
            plt.plot( Xk1, Yk1,                       # 画出第 1 层可移动格点 | Plot the level-1 moving grid points
                      color='blue',                   # 全部点设置为蓝色 | Set all points to blue
                      marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
    
    if (input_data.moving_grid_level > 2):
        for k in range(input_data.puncture_number):
            Xk2, Yk2 = numpy.meshgrid( initial_AMR_grid_data.Moving_Grid_X[2,k], initial_AMR_grid_data.Moving_Grid_Y[2,k] )
            plt.plot( Xk2, Yk2,                       # 画出第 1 层可移动格点 | Plot the level-1 moving grid points
                      color='navy',                   # 全部点设置为深蓝色 | Set all points to dark blue
                      marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines

    if (input_data.moving_grid_level > 3):
        for k in range(input_data.puncture_number):
            Xk3, Yk3 = numpy.meshgrid( initial_AMR_grid_data.Moving_Grid_X[3,k], initial_AMR_grid_data.Moving_Grid_Y[3,k] )
            plt.plot( Xk3, Yk3,                       # 画出第 1 层可移动格点 | Plot the level-1 moving grid points
                      color='gray',                   # 全部点设置为灰色 | Set all points to gray
                      marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
    
    if (input_data.moving_grid_level > 4):
        for k in range(input_data.puncture_number):
            Xk4, Yk4 = numpy.meshgrid( initial_AMR_grid_data.Moving_Grid_X[4,k], initial_AMR_grid_data.Moving_Grid_Y[4,k] )
            plt.plot( Xk4, Yk4,                       # 画出第 1 层可移动格点 | Plot the level-1 moving grid points
                      color='black',                  # 全部点设置为灰色 | Set all points to gray
                      marker='.',                     # 点的形状为圆点 | The marker shape is a round dot
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接 | Empty line style, i.e. the points are not connected by lines
    
    plt.grid(True)
    ## plt.show()
    plt.savefig( os.path.join(input_data.File_directory, "Initial_Grid.jpeg") )
    plt.savefig( os.path.join(input_data.File_directory, "Initial_Grid.pdf")  )
    
    return

#########################################################################



#########################################################################
    
## 该函数将以上格点信息写入 AMSS-NCKU 程序的输入文件 | This function writes the above grid information into the input file of the AMSS-NCKU program

#########################################################################
    
def append_AMSSNCKU_cgh_input( input_language, initial_AMR_grid_data ): 

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU.input"), "a")  
    # "a" 表示追加输出 | "a" means appending to the output

    ## 输出 cgh 的相关设定 | Output the cgh-related settings

    print( file=file1 )
    print( "cgh::moving levels start from = ", input_data.static_grid_level, file=file1 )
    print( "cgh::levels = ",                   input_data.grid_level,        file=file1 )

    ## 输出固定网格的信息 | Output the static grid information

    for i in range(input_data.static_grid_level): 

        print( f"cgh::grids[{i}]       = 1",                                file=file1 )
        
        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::shape[{i}][0][0] = { initial_AMR_grid_data.static_grid_number_x//2 } ", file=file1 )
            print( f"cgh::shape[{i}][0][1] = { initial_AMR_grid_data.static_grid_number_y//2 } ", file=file1 )
        else:
            print( f"cgh::shape[{i}][0][0] = { initial_AMR_grid_data.static_grid_number_x } ",    file=file1 )
            print( f"cgh::shape[{i}][0][1] = { initial_AMR_grid_data.static_grid_number_y } ",    file=file1 )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { initial_AMR_grid_data.static_grid_number_z//2 } ", file=file1 )
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { initial_AMR_grid_data.static_grid_number_z//2 } ", file=file1 )
        elif ( input_data.Symmetry == "no-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { initial_AMR_grid_data.static_grid_number_z } ",    file=file1 )
        else:
            print( " Symmetry Setting Error " )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::bbox[{i}][0][0]  = 0.0 ", file=file1 )
            print( f"cgh::bbox[{i}][0][1]  = 0.0 ", file=file1 )
        else:
            print( f"cgh::bbox[{i}][0][0]  = { initial_AMR_grid_data.Grid_X_Min[i] } ", file=file1 )
            print( f"cgh::bbox[{i}][0][1]  = { initial_AMR_grid_data.Grid_Y_Min[i] } ", file=file1 )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = 0.0 ", file=file1 )
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = 0.0 ", file=file1 )
        elif ( input_data.Symmetry == "no-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = { initial_AMR_grid_data.Grid_Z_Min[i] } ", file=file1 )
        else:
            print( " Symmetry Setting Error " )

        print( f"cgh::bbox[{i}][0][3]  = { initial_AMR_grid_data.Grid_X_Max[i] } ", file=file1 )
        print( f"cgh::bbox[{i}][0][4]  = { initial_AMR_grid_data.Grid_Y_Max[i] } ", file=file1 )
        print( f"cgh::bbox[{i}][0][5]  = { initial_AMR_grid_data.Grid_Z_Max[i] } ", file=file1 )

    ## 输出可移动网格的信息 | Output the moving grid information 
       
    ## 对可移动网格层数进行循环 | Loop over the moving grid levels
    
    for i in range(input_data.moving_grid_level):

        j = i + input_data.static_grid_level
        print( f"cgh::grids[{j}]       = { input_data.puncture_number }", file=file1 )

        ## 对 Puncture 数目进行循环 | Loop over the number of Punctures
        for k in range(input_data.puncture_number):

            if ( input_data.Symmetry == "octant-symmetry" ): 
                print( f"cgh::shape[{j}][{k}][0] = { input_data.moving_grid_number//2 } ", file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { input_data.moving_grid_number//2 } ", file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { input_data.moving_grid_number//2 } ", file=file1 )
            elif ( input_data.Symmetry == "equatorial-symmetry" ): 
                print( f"cgh::shape[{j}][{k}][0] = { input_data.moving_grid_number } ",    file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { input_data.moving_grid_number } ",    file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { input_data.moving_grid_number//2 } ", file=file1 )
            elif ( input_data.Symmetry == "no-symmetry" ):
                print( f"cgh::shape[{j}][{k}][0] = { input_data.moving_grid_number } ",    file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { input_data.moving_grid_number } ",    file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { input_data.moving_grid_number } ",    file=file1 )   
            else:
                print( " Symmetry Setting Error" )

            ## 原设定 | Original setup
            ## print( f"cgh::bbox[{j}][{k}][0]  = { initial_AMR_grid_data.Moving_Grid_X_Min[i,k] } ", file=file1 )
            ## print( f"cgh::bbox[{j}][{k}][1]  = { initial_AMR_grid_data.Moving_Grid_Y_Min[i,k] } ", file=file1 )

            if ( input_data.Symmetry == "octant-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { max(0.0, initial_AMR_grid_data.Moving_Grid_X_Min[i,k]) } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { max(0.0, initial_AMR_grid_data.Moving_Grid_Y_Min[i,k]) } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { max(0.0, initial_AMR_grid_data.Moving_Grid_Z_Min[i,k]) } ", file=file1 )
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { initial_AMR_grid_data.Moving_Grid_X_Min[i,k] } ",           file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { initial_AMR_grid_data.Moving_Grid_Y_Min[i,k] } ",           file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { max(0.0, initial_AMR_grid_data.Moving_Grid_Z_Min[i,k]) } ", file=file1 )
            elif ( input_data.Symmetry == "no-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { initial_AMR_grid_data.Moving_Grid_X_Min[i,k] } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { initial_AMR_grid_data.Moving_Grid_Y_Min[i,k] } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { initial_AMR_grid_data.Moving_Grid_Z_Min[i,k] } ", file=file1 )
            else:
                print( " Symmetry Setting Error" )

            print( f"cgh::bbox[{j}][{k}][3]  = { initial_AMR_grid_data.Moving_Grid_X_Max[i,k] } ",     file=file1 )
            print( f"cgh::bbox[{j}][{k}][4]  = { initial_AMR_grid_data.Moving_Grid_Y_Max[i,k] } ",     file=file1 )
            print( f"cgh::bbox[{j}][{k}][5]  = { initial_AMR_grid_data.Moving_Grid_Z_Max[i,k] } ",     file=file1 )

    ## 输出 BSSN 的相关设定 | Output the BSSN-related settings
    
    print(                                                                         file=file1 )
    print( "############ for shell-box coupling set this exactly to box boundary", file=file1 )
    print( "BSSN::Shell shape[0]   = ", input_data.shell_grid_number[0],           file=file1 )
    print( "BSSN::Shell shape[1]   = ", input_data.shell_grid_number[1],           file=file1 )
    print( "BSSN::Shell shape[2]   = ", input_data.shell_grid_number[2],           file=file1 )
    print( "BSSN::Shell R range[0] = ", initial_AMR_grid_data.Shell_R_Min,         file=file1 )
    print( "BSSN::Shell R range[1] = ", initial_AMR_grid_data.Shell_R_Max,         file=file1 )
    print(                                                                         file=file1 )
            
    file1.close()
    
    ## 屏幕输出： cgh 的信息已写入 AMSS-NCKU 输入文件 | Screen output: the cgh information has been written into the AMSS-NCKU input file
    if ( input_language == "Chinese" ):
        print(                                                               )
        print( " 完成 AMSS-NCKU 程序的输入文件 "                                )    
        print( " 但是 AMSS-NCKU ABE 中 TwoPuncture 相关的输入要后几步运行后追加 "  )
        print(                                                               )
    elif ( input_language == "English" ): 
        print(                                                                                 )
        print( " The input parfile for AMSS-NCKU C++ executable file ABE has been generated. " )    
        print( " However, the input relevant to TwoPuncture need to be appended later. "       )
        print(                                                                                 )

    return 
    
#########################################################################


    
