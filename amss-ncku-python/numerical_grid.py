
#################################################
##
## 这个文件包含了数值相对论所需要的格点
## 小曲
## 2024/03/20
## 2025/09/14 修改 
##
#################################################

import numpy                              ## 导入 numpy 包 
import matplotlib.pyplot as plt           ## 导入 matplotlib 包中的 pyplot 函数，并用 plt 名称代表 matplotlib.pyplot
import os                                 ## 导入 os 包进行系统操作

import AMSS_NCKU_Input   as input_data    ## 导入程序输入文件
## import print_information

#################################################

# 设置黑洞 puncture （穿刺法）的信息

puncture = numpy.zeros( (input_data.puncture_number,3) )      # 初始化每个 puncture 的位置

print(                                   )
print( " 正在设定 Puncture 位置和动量参数" )
print( " Setting Puncture's position and momentum " )
print(                                              )

#################################################

## 设置 puncture 的位置

## 如果设置 puncture 轨道位置的方式选为 Automatically-BBH，则读入重设后的 puncture 位置 
 
if (input_data.puncture_data_set == "Automatically-BBH" ):

    import generate_TwoPuncture_input
    
    for i in range(input_data.puncture_number):
        if (i<=1):
            puncture[i] = generate_TwoPuncture_input.position_BH[i]
        else:
            puncture[i] = input_data.position_BH[i]
    
## 如果如果设置 puncture 轨道位置的方式选为 Manually，则直接读入  puncture 位置  
   
elif (input_data.puncture_data_set == "Manually" ):

    puncture = input_data.position_BH 
    
else: 
   
   print(                                          )
   print( " 正在设定 Puncture 位置和动量参数的设定错误" )
   print( " Found Error in setting Puncture's position and momentum !!! " )
   print(                                                                 )

#################################################

## 输出网格的信息

print(                                     )   
print( " Wirte Down The Grid Information " )
print(                                     )   
print( " Number of Total Grid Level = ",          input_data.grid_level        )      ## 输出网格的总层数
print( " Number of Static Grid Level = ",         input_data.static_grid_level )      ## 输出静态网格的层数
print( " Number of Moving Grid Level = ",         input_data.moving_grid_level )      ## 输出可移动网格的层数
## print( " Number of Points in Each Grid Level = ", input_data.grid_number    )      ## 输出每一层每一维数的网格格点数目
print(                                                                         )

#################################################

print(                     )
print( " 正在设定计算网格 " )
print( " Setting the demanded numerical grid " )
print(                                         )

#################################################

## 初始化网格信息

## 初始化网格坐标最小值、最大值、格点数目

Grid_X_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 X 坐标的最小值
Grid_X_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 X 坐标的最大值
Grid_Y_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Y 坐标的最小值
Grid_Y_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Y 坐标的最大值
Grid_Z_Min = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Z 坐标的最小值
Grid_Z_Max = numpy.zeros( (input_data.grid_level) )    # 定义实数数组，作为每一层网格 Z 坐标的最大值

Grid_Resolution = numpy.zeros( input_data.grid_level ) # 定义实数数组，作为每层网格分辨率

largest_box_X_Max = input_data.largest_box_xyz_max[0] 
largest_box_Y_Max = input_data.largest_box_xyz_max[1]
largest_box_Z_Max = input_data.largest_box_xyz_max[2]
largest_box_X_Min = input_data.largest_box_xyz_min[0] 
largest_box_Y_Min = input_data.largest_box_xyz_min[1]
largest_box_Z_Min = input_data.largest_box_xyz_min[2]

# 定义整数数组，作为每一层固定网格每一维的格点数
static_grid_number_x = input_data.static_grid_number 
static_grid_number_y = int( (largest_box_Y_Max - largest_box_Y_Min) * ( static_grid_number_x / (largest_box_X_Max-largest_box_X_Min) ) )
static_grid_number_z = int( (largest_box_Z_Max - largest_box_Z_Min) * ( static_grid_number_x / (largest_box_X_Max-largest_box_X_Min) ) )
    
# 定义整数数组，作为每一层可移动网格每一维的格点数
moving_grid_number   = input_data.moving_grid_number   

#################################################

## 初始化固定网格

# 将每个方向的网格数目自动调整为偶数
# Python 中使用 % 取余数
# print(static_grid_number_x % 2)
if ( (static_grid_number_x % 2) != 0) :
    static_grid_number_x = static_grid_number_x + 1
if ( (static_grid_number_y % 2) != 0) :
    static_grid_number_y = static_grid_number_y + 1
if ( (static_grid_number_z % 2) != 0) :
    static_grid_number_z = static_grid_number_z + 1
# 为了防止移动网格移动时边界无法与固定网格对其，我们进一步要求每个方向的网格数目自动调整为 4 的倍数
if ( (static_grid_number_x % 4) != 0) :
    static_grid_number_x = static_grid_number_x + 2
if ( (static_grid_number_y % 4) != 0) :
    static_grid_number_y = static_grid_number_y + 2
if ( (static_grid_number_z % 4) != 0) :
    static_grid_number_z = static_grid_number_z + 2
'''
# 为了防止移动网格移动时边界无法与固定网格对其，我们进一步要求每个方向的网格数目自动调整为 8 的倍数
if ( (static_grid_number_x % 8) != 0) :
    static_grid_number_x = static_grid_number_x + 4
if ( (static_grid_number_y % 8) != 0) :
    static_grid_number_y = static_grid_number_y + 4
if ( (static_grid_number_z % 8) != 0) :
    static_grid_number_z = static_grid_number_z + 4
'''

##  定义实数数组，维度 grid_number * static_grid_level，分别作为每一层固定网格的 X Y Z 坐标
Static_Grid_X = numpy.zeros( (input_data.static_grid_level, static_grid_number_x+1) )   
Static_Grid_Y = numpy.zeros( (input_data.static_grid_level, static_grid_number_y+1) )   
Static_Grid_Z = numpy.zeros( (input_data.static_grid_level, static_grid_number_z+1) )  

#################################################

## 初始化可移动网格

##  定义实数数组，维度 grid_number * puncture_number * moving_grid_level，作为每一层可移动网格的 X Y Z 坐标
Moving_Grid_X = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) )
Moving_Grid_Y = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) ) 
Moving_Grid_Z = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number, input_data.moving_grid_number+1) )

#################################################

## 初始化可移动网格的坐标最小值、最大值

Moving_Grid_X_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) ) 
Moving_Grid_X_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
Moving_Grid_Y_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
Moving_Grid_Y_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )                
Moving_Grid_Z_Min = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )
Moving_Grid_Z_Max = numpy.zeros( (input_data.moving_grid_level, input_data.puncture_number) )          

#################################################

## 设置每层网格的网格分辨率

for i in range(input_data.static_grid_level) :
    if i==0:
        Grid_Resolution[i] = ( largest_box_X_Max - largest_box_X_Min ) / static_grid_number_x
    else:
        Grid_Resolution[i] = Grid_Resolution[i-1] / input_data.devide_factor
    
for j in range(input_data.moving_grid_level) : 
    i = j + input_data.static_grid_level
    Grid_Resolution[i] = Grid_Resolution[i-1] / input_data.devide_factor

#################################################

## 根据输入文件，设置每层固定网格坐标的最小值和最大值

## 设置第一层固定网格的坐标最小值与最大值
Grid_X_Min[0] = largest_box_X_Min 
Grid_X_Max[0] = largest_box_X_Max
Grid_Y_Min[0] = largest_box_Y_Min  
Grid_Y_Max[0] = largest_box_Y_Max
Grid_Z_Min[0] = largest_box_Z_Min
Grid_Z_Max[0] = largest_box_Z_Max
## 重新调整 yx 的最大坐标，确保 xyz 三个方向的分辨率一致
Grid_Y_Max[0] = Grid_Y_Min[0] + Grid_Resolution[0] * static_grid_number_y
Grid_Z_Max[0] = Grid_Z_Min[0] + Grid_Resolution[0] * static_grid_number_z
## 如果系统存在对称性，还要进行进一步的调整
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

print( " 正在调整固定网格格点的位置，使得坐标原点位于固定网格上 " )
print( " adjusting the static gird points, making the original point (0,0,0) to the static gird points " )
print(                                                                                                   )

## 设置其它层固定网格的坐标最大值与最小值
for i in range(input_data.static_grid_level-1) :
    ## 如果坐标原点不在最外层固定网格上，调整最外层固定网格，使得原点位于最外层网格上
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
    ## 每层固定网格最小值（或最大值）为上一层固定网格最小值（或最大值）除以因子 devide_factor
    Grid_X_Min[i+1] = Grid_X_Min[i] / input_data.devide_factor    
    Grid_X_Max[i+1] = Grid_X_Max[i] / input_data.devide_factor    
    Grid_Y_Min[i+1] = Grid_Y_Min[i] / input_data.devide_factor
    Grid_Y_Max[i+1] = Grid_Y_Max[i] / input_data.devide_factor
    Grid_Z_Min[i+1] = Grid_Z_Min[i] / input_data.devide_factor
    Grid_Z_Max[i+1] = Grid_Z_Max[i] / input_data.devide_factor
    

## 固定网格与可移动网格的格点数目不同，为保证相同分辨率，需要引入网格大小调节因子
adjust_factor = input_data.moving_grid_number / input_data.static_grid_number

## 设置第一层可移动网格的坐标最大值与最小值
i = input_data.static_grid_level 
if (i < input_data.grid_level):
    Grid_X_Min[i] = ( Grid_X_Min[i-1] / input_data.devide_factor ) * adjust_factor
    Grid_X_Max[i] = - Grid_X_Min[i]
    # Grid_X_Max[i] = ( Grid_X_Max[i-1] / input_data.devide_factor ) * adjust_factor   
    ## 原来的设定
    # Grid_Y_Min[i] = ( Grid_Y_Min[i-1] / input_data.devide_factor ) * adjust_factor
    # Grid_Y_Max[i] = ( Grid_Y_Max[i-1] / input_data.devide_factor ) * adjust_factor
    # Grid_Z_Min[i] = ( Grid_Z_Min[i-1] / input_data.devide_factor ) * adjust_factor
    # Grid_Z_Max[i] = ( Grid_Z_Max[i-1] / input_data.devide_factor ) * adjust_factor
    ## 现在的设定，总是保证移动网格为立方网格
    Grid_Y_Min[i] = Grid_X_Min[i]
    Grid_Y_Max[i] = Grid_X_Max[i]
    Grid_Z_Min[i] = Grid_X_Min[i]
    Grid_Z_Max[i] = Grid_X_Max[i]

    # print( " Grid_X_Max[i] = ", Grid_X_Max[i] )
    # print( " Grid_Y_Max[i] = ", Grid_Y_Max[i] )

## 设置其它层可移动网格的坐标最大值与最小值
for j in range(input_data.moving_grid_level-1) :
    k = input_data.static_grid_level + j
    Grid_X_Min[k+1] = Grid_X_Min[k] / input_data.devide_factor    
    Grid_X_Max[k+1] = Grid_X_Max[k] / input_data.devide_factor    
    Grid_Y_Min[k+1] = Grid_Y_Min[k] / input_data.devide_factor
    Grid_Y_Max[k+1] = Grid_Y_Max[k] / input_data.devide_factor
    Grid_Z_Min[k+1] = Grid_Z_Min[k] / input_data.devide_factor
    Grid_Z_Max[k+1] = Grid_Z_Max[k] / input_data.devide_factor
    
## 设置最外层球壳格点 shell patch 的坐标最大值和最小值

Shell_R_Resolution = Grid_Resolution[0]
Shell_R_Min        = largest_box_X_Max
Shell_R_Max        = Shell_R_Min + Grid_Resolution[0] * input_data.shell_grid_number[2]

#################################################


#################################################

## 设置每层网格的坐标值

#################################################

## 设置固定网格的格点坐标坐标

## 线性格点

if input_data.static_grid_type == 'Linear' : 

    for i in range(input_data.static_grid_level):
        Static_Grid_X[i] = numpy.linspace( Grid_X_Min[i], Grid_X_Max[i], static_grid_number_x+1 ) 
        Static_Grid_Y[i] = numpy.linspace( Grid_Y_Min[i], Grid_Y_Max[i], static_grid_number_y+1 ) 
        Static_Grid_Z[i] = numpy.linspace( Grid_Z_Min[i], Grid_Z_Max[i], static_grid_number_z+1 )
     # 利用 numpy 设置线性格点，参数表示起始值 Rmin，最大值 Rmax，格点数目 Rnum
     # 注意，如果是线性格点，则格点坐标最大值为 GridMax；如果是对数格点，则格点坐标最大值为 e^{GridMax}

else:
    print(                                                       )
    print( " Static Grid Error: Grid Type is Undifined !!!!!!! " )
    print(                                                       )

#################################################
  
## 设置可移动网格的格点坐标坐标  

print(                                                                                                            )
print( " 正在调整移动网格中心的位置，以便同时与 Puncture 位置以及上一层粗网格相匹配 "                                       )
print( " adjusting the moving gird points, ensuring the alliance of moving grids points and static grids points " )
print(                                                                                                            )

## 为了让每一层可移动网格的边界与上一层网格重合，要对 Puncture 的位置进行调整
adjust_puncture = numpy.zeros( (input_data.puncture_number, 3) )

## 线性格点
  
if ( input_data.moving_grid_type == "Linear" ): 
    
    ## 对可移动网格层数进行循环
    for j in range(input_data.moving_grid_level) : 

        i = j + input_data.static_grid_level
        
        ## 对 Puncture 数目进行循环
        for k in range(input_data.puncture_number) :
            
            ## 将 Puncture 的位置进行调整
            if j==0 :
            
                level0 = input_data.static_grid_level - 1  ## 引入新变量，避免代码过长
                
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
            ## Puncture 的位置调整完毕
            
            ## 为了防止 C++ 程度读入数据出错
            ## 如果遇到很小的数 1e-10，保留 2 位小数，不使用科学技术法
            if ( abs(adjust_puncture[k,0]) < 1e-10 ):
                adjust_puncture[k,0] = 0.00
                # adjust_puncture[k,0] = f"{ adjust_puncture[k,0]:.2f }
            if ( abs(adjust_puncture[k,1]) < 1e-10 ):
                adjust_puncture[k,1] = 0.00
                # adjust_puncture[k,1] = f"{ adjust_puncture[k,1]:.2f }
            if ( abs(adjust_puncture[k,2]) < 1e-10 ):
                adjust_puncture[k,2] = 0.00
                # adjust_puncture[k,2] = f"{ adjust_puncture[k,2]:.2f }
            
            # 第 j 层可移动网格的 XYZ 坐标最小值（或最大值） = 第 i 层网格的 XYZ 坐标最小值（或最大值） + 第 k 个 puncture 的 XYZ 坐标
            Moving_Grid_X_Min[j,k] = adjust_puncture[k,0] + Grid_X_Min[i]  
            Moving_Grid_X_Max[j,k] = adjust_puncture[k,0] + Grid_X_Max[i]  
            Moving_Grid_Y_Min[j,k] = adjust_puncture[k,1] + Grid_Y_Min[i]  
            Moving_Grid_Y_Max[j,k] = adjust_puncture[k,1] + Grid_Y_Max[i]  
            Moving_Grid_Z_Min[j,k] = adjust_puncture[k,2] + Grid_Z_Min[i]  
            Moving_Grid_Z_Max[j,k] = adjust_puncture[k,2] + Grid_Z_Max[i]
            
            ## 为了防止 C++ 程度读入数据出错
            ## 如果遇到很小的数 1e-10，保留 2 位小数，不使用科学技术法
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
            
            print( f" adjust_puncture[{i},{k},0] = { adjust_puncture[k,0] } " )
            print( f" adjust_puncture[{i},{k},1] = { adjust_puncture[k,1] } " )
            print( f" adjust_puncture[{i},{k},2] = { adjust_puncture[k,2] } " )
            
            ## 利用 numpy 设置线性格点，参数表示起始值 Rmin，最大值 Rmax，格点数目 Rnum
            Moving_Grid_X[j,k] = numpy.linspace( Moving_Grid_X_Min[j,k], Moving_Grid_X_Max[j,k], moving_grid_number + 1 )  
            Moving_Grid_Y[j,k] = numpy.linspace( Moving_Grid_Y_Min[j,k], Moving_Grid_Y_Max[j,k], moving_grid_number + 1 )
            Moving_Grid_Z[j,k] = numpy.linspace( Moving_Grid_Z_Min[j,k], Moving_Grid_Z_Max[j,k], moving_grid_number + 1 )
        
else:
    print(                                                       )
    print( " Moving Grid Error: Grid Type is Undifined !!!!!!! " )
    print(                                                       )

print(                            )
print( " 移动网格中心位置调整完毕 " )
print(                            )

#################################################


#################################################

## 该函数画出初始的网格
    
def plot_initial_grid():

## 依次画出最后的网格

    if (input_data.static_grid_level > 0):
        X0, Y0 = numpy.meshgrid( Static_Grid_X[0], Static_Grid_Y[0] )
        plt.plot( X0, Y0,                         # 画出第 0 层静态格点
                  color='brown',  	          # 全部点设置为红棕色
                  marker='.',  	                  # 点的形状为圆点
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接
              
    if (input_data.static_grid_level > 1):
        X1, Y1 = numpy.meshgrid( Static_Grid_X[1], Static_Grid_Y[1] )
        plt.plot( X1, Y1,                         # 画出第 1 层静态格点
                  color='red',                    # 全部点设置为红色
                  marker='.',                     # 点的形状为圆点
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接
                  
    if (input_data.static_grid_level > 2):
        X2, Y2 = numpy.meshgrid( Static_Grid_X[2], Static_Grid_Y[2] )
        plt.plot( X2, Y2,                         # 画出第 2 层静态格点
                  color='orange',                 # 全部点设置为橙色
                  marker='.',                     # 点的形状为圆点
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接
                  
    if (input_data.static_grid_level > 3):
        X3, Y3 = numpy.meshgrid( Static_Grid_X[3], Static_Grid_Y[3] )
        plt.plot( X3, Y3,                         # 画出第 3 层静态格点
                  color='yellow',                 # 全部点设置为黄色
                  marker='.',                     # 点的形状为圆点
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接
                  
    if (input_data.static_grid_level > 4):
        X4, Y4 = numpy.meshgrid( Static_Grid_X[4], Static_Grid_Y[4] )
        plt.plot( X4, Y4,                         # 画出第 4 层静态格点
                  color='greenyellow',            # 全部点设置为黄绿色
                  marker='.',                     # 点的形状为圆点
                  linestyle='' )                  # 线型为空，也即点与点之间不用线连接
                  
    ## 依次画出可移动网格

    if (input_data.moving_grid_level > 0):
        for k in range(input_data.puncture_number):
            Xk0, Yk0 = numpy.meshgrid( Moving_Grid_X[0,k], Moving_Grid_Y[0,k] )
            plt.plot( Xk0, Yk0,                       # 画出第 0 层可移动格点
                      color='cyan',  	              # 全部点设置为青蓝色
                      marker='.',  	              # 点的形状为圆点
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接
 
    if (input_data.moving_grid_level > 1):
        for k in range(input_data.puncture_number):
            Xk1, Yk1 = numpy.meshgrid( Moving_Grid_X[1,k], Moving_Grid_Y[1,k] )
            plt.plot( Xk1, Yk1,                       # 画出第 1 层可移动格点
                      color='blue',                   # 全部点设置为蓝色
                      marker='.',                     # 点的形状为圆点
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接
    
    if (input_data.moving_grid_level > 2):
        for k in range(input_data.puncture_number):
            Xk2, Yk2 = numpy.meshgrid( Moving_Grid_X[2,k], Moving_Grid_Y[2,k] )
            plt.plot( Xk2, Yk2,                       # 画出第 1 层可移动格点
                      color='navy',                   # 全部点设置为深蓝色
                      marker='.',                     # 点的形状为圆点
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接

    if (input_data.moving_grid_level > 3):
        for k in range(input_data.puncture_number):
            Xk3, Yk3 = numpy.meshgrid( Moving_Grid_X[3,k], Moving_Grid_Y[3,k] )
            plt.plot( Xk3, Yk3,                       # 画出第 1 层可移动格点
                      color='gray',                   # 全部点设置为灰色
                      marker='.',                     # 点的形状为圆点
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接
    
    if (input_data.moving_grid_level > 4):
        for k in range(input_data.puncture_number):
            Xk4, Yk4 = numpy.meshgrid( Moving_Grid_X[4,k], Moving_Grid_Y[4,k] )
            plt.plot( Xk4, Yk4,                       # 画出第 1 层可移动格点
                      color='black',                  # 全部点设置为灰色
                      marker='.',                     # 点的形状为圆点
                      linestyle='' )                  # 线型为空，也即点与点之间不用线连接
    
    plt.grid(True)
    ## plt.show()
    plt.savefig( os.path.join(input_data.File_directory, "Initial_Grid.jpeg") )
    plt.savefig( os.path.join(input_data.File_directory, "Initial_Grid.pdf")  )

#################################################


#################################################
    
## 将以上格点信息写入 AMSS-NCKU 程序的输入文件
    
def append_AMSSNCKU_cgh_input(): 

    file1 = open( os.path.join(input_data.File_directory, "AMSS-NCKU.input"), "a")  
    # "a" 表示追加输出

    ## 输出 cgh 的相关设定

    print( file=file1 )
    print( "cgh::moving levels start from = ", input_data.static_grid_level, file=file1 )
    print( "cgh::levels = ",                   input_data.grid_level,        file=file1)

    ## 输出固定网格的信息

    for i in range(input_data.static_grid_level): 

        print( f"cgh::grids[{i}]       = 1",                                file=file1 )
        
        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::shape[{i}][0][0] = { static_grid_number_x//2 } ", file=file1 )
            print( f"cgh::shape[{i}][0][1] = { static_grid_number_y//2 } ", file=file1 )
        else:
            print( f"cgh::shape[{i}][0][0] = { static_grid_number_x } ",    file=file1 )
            print( f"cgh::shape[{i}][0][1] = { static_grid_number_y } ",    file=file1 )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { static_grid_number_z//2 } ", file=file1 )
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { static_grid_number_z//2 } ", file=file1 )
        elif ( input_data.Symmetry == "no-symmetry" ):
            print( f"cgh::shape[{i}][0][2] = { static_grid_number_z } ",    file=file1 )
        else:
            print( " Symmetry Setting Error " )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::bbox[{i}][0][0]  = 0.0 ",                       file=file1 )
            print( f"cgh::bbox[{i}][0][1]  = 0.0 ",                       file=file1 )
        else:
            print( f"cgh::bbox[{i}][0][0]  = { Grid_X_Min[i] } ",         file=file1 )
            print( f"cgh::bbox[{i}][0][1]  = { Grid_Y_Min[i] } ",         file=file1 )

        if ( input_data.Symmetry == "octant-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = 0.0 ",                       file=file1 )
        elif ( input_data.Symmetry == "equatorial-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = 0.0 ",                       file=file1 )
        elif ( input_data.Symmetry == "no-symmetry" ):
            print( f"cgh::bbox[{i}][0][2]  = { Grid_Z_Min[i] } ",         file=file1 )
        else:
            print( " Symmetry Setting Error " )

        print( f"cgh::bbox[{i}][0][3]  = { Grid_X_Max[i] } ",             file=file1 )
        print( f"cgh::bbox[{i}][0][4]  = { Grid_Y_Max[i] } ",             file=file1 )
        print( f"cgh::bbox[{i}][0][5]  = { Grid_Z_Max[i] } ",             file=file1 )

    ## 输出可移动网格的信息 
       
    ## 对可移动网格层数进行循环
    
    for i in range(input_data.moving_grid_level):

        j = i + input_data.static_grid_level
        print( f"cgh::grids[{j}]       = { input_data.puncture_number }",        file=file1 )

        ## 对 Puncture 数目进行循环
        for k in range(input_data.puncture_number):

            if ( input_data.Symmetry == "octant-symmetry" ): 
                print( f"cgh::shape[{j}][{k}][0] = { moving_grid_number//2 } ",  file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { moving_grid_number//2 } ",  file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { moving_grid_number//2 } ",  file=file1 )
            elif ( input_data.Symmetry == "equatorial-symmetry" ): 
                print( f"cgh::shape[{j}][{k}][0] = { moving_grid_number } ",     file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { moving_grid_number } ",     file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { moving_grid_number//2 } ",  file=file1 )
            elif ( input_data.Symmetry == "no-symmetry" ):
                print( f"cgh::shape[{j}][{k}][0] = { moving_grid_number } ",     file=file1 )
                print( f"cgh::shape[{j}][{k}][1] = { moving_grid_number } ",     file=file1 )
                print( f"cgh::shape[{j}][{k}][2] = { moving_grid_number } ",     file=file1 )   
            else:
                print( " Symmetry Setting Error" )

            print( f"cgh::bbox[{j}][{k}][0]  = { Moving_Grid_X_Min[i,k] } ",     file=file1 )
            print( f"cgh::bbox[{j}][{k}][1]  = { Moving_Grid_Y_Min[i,k] } ",     file=file1 )

            if ( input_data.Symmetry == "octant-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { max(0.0, Moving_Grid_X_Min[i,k]) } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { max(0.0, Moving_Grid_Y_Min[i,k]) } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { max(0.0, Moving_Grid_Z_Min[i,k]) } ", file=file1 )
            elif ( input_data.Symmetry == "equatorial-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { Moving_Grid_X_Min[i,k] } ",           file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { Moving_Grid_Y_Min[i,k] } ",           file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { max(0.0, Moving_Grid_Z_Min[i,k]) } ", file=file1 )
            elif ( input_data.Symmetry == "no-symmetry" ):
                print( f"cgh::bbox[{j}][{k}][0]  = { Moving_Grid_X_Min[i,k] } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][1]  = { Moving_Grid_Y_Min[i,k] } ", file=file1 )
                print( f"cgh::bbox[{j}][{k}][2]  = { Moving_Grid_Z_Min[i,k] } ", file=file1 )
            else:
                print( " Symmetry Setting Error" )

            print( f"cgh::bbox[{j}][{k}][3]  = { Moving_Grid_X_Max[i,k] } ",     file=file1 )
            print( f"cgh::bbox[{j}][{k}][4]  = { Moving_Grid_Y_Max[i,k] } ",     file=file1 )
            print( f"cgh::bbox[{j}][{k}][5]  = { Moving_Grid_Z_Max[i,k] } ",     file=file1 )

    ## 输出 BSSN 的相关设定
    
    print(                                                                         file=file1 )
    print( "############ for shell-box coupling set this exactly to box boundary", file=file1 )
    print( "BSSN::Shell shape[0]   = ", input_data.shell_grid_number[0],           file=file1 )
    print( "BSSN::Shell shape[1]   = ", input_data.shell_grid_number[1],           file=file1 )
    print( "BSSN::Shell shape[2]   = ", input_data.shell_grid_number[2],           file=file1 )
    print( "BSSN::Shell R range[0] = ", Shell_R_Min,                               file=file1 )
    print( "BSSN::Shell R range[1] = ", Shell_R_Max,                               file=file1 )
    print(                                                                         file=file1 )
            
    file1.close()

    return file1
    
#################################################


    
