
#################################################
##
## 这个文件对 AMSS-NCKU 数值相对论的结果进行画图
## 小曲
## 2024/10/01 --- 2025/09/14
##
#################################################

import numpy                               ## 导入 numpy 包进行数组的操作
import matplotlib.pyplot    as     plt     ## 导入 matplotlib 包进行画图
from   mpl_toolkits.mplot3d import Axes3D  ## 画 3 维图需要
import glob
import os                                  ## 导入 os 包进行系统操作

import plot_binary_data
import AMSS_NCKU_Input as input_data

# plt.rcParams['text.usetex'] = True  ## 在绘图中允许使用 latex 字体



####################################################################################

## 该函数根据二进制数据画出所有二维图

def generate_binary_data_plot( binary_outdir, figure_outdir ):

    # 生成若干文件夹存放图片

    surface_plot_outdir = os.path.join( figure_outdir, "surface plot" )
    os.mkdir( surface_plot_outdir )

    density_plot_outdir = os.path.join( figure_outdir, "density plot" )
    os.mkdir( density_plot_outdir )

    contour_plot_outdir = os.path.join( figure_outdir, "contour plot" )
    os.mkdir( contour_plot_outdir )

    print(                                  )
    print( " 读取 AMSS-NCKU 程序的二进制数据 " )
    print( " Reading AMSS-NCKU Binary Data From Output " )
    print(                                               )

    print( " 二进制数据列表 " )
    
    ## 设置对什么文件画图（这里设置对所有二进制文件画图）
    globby = glob.glob( os.path.join(binary_outdir, '*.bin') ) 
    file_list = []
    for x in sorted(globby):
        file_list.append(x)
        print(x)

    ## 对列表中的所有文件画图
    for filename in file_list:
        print(filename)
        plot_binary_data.plot_binary_data(filename, binary_outdir, figure_outdir)

    print(                        )
    print( " 二进制数据画图已完成 " )
    print( " Binary Data Plot Has been Finished " )
    print(                                        )

    return

####################################################################################



####################################################################################

## 该函数对黑洞轨迹画图

def generate_puncture_orbit_plot( outdir, figure_outdir ):

    print(                                                    )
    print( " 正在对黑洞轨迹进行画图 "                            )
    print( " Plotting the black holes' trajectory (2D plot)" )
    print(                                                   )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    print( " 对应数据文件为 ",              file0 )
    print( " Corresponding data file = ", file0 )

    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)

    # print(data[:,0])
    # print(data[:,2])

    # 初始化黑洞坐标的最大值和最小值
    BH_Xmin = numpy.zeros(input_data.puncture_number)
    BH_Xmax = numpy.zeros(input_data.puncture_number)
    BH_Ymin = numpy.zeros(input_data.puncture_number)
    BH_Ymax = numpy.zeros(input_data.puncture_number)
    BH_Zmin = numpy.zeros(input_data.puncture_number)
    BH_Zmax = numpy.zeros(input_data.puncture_number)
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（XY图）
    
    plt.figure( figsize=(8,8)                         )   ## 这里 figsize 可以设定图形的大小
    plt.title( " Black Hole Trajectory ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    
    for i in range(input_data.puncture_number):
        BH_x       = data[:, 3*i+1]
        BH_y       = data[:, 3*i+2]
        BH_z       = data[:, 3*i+3]
        BH_Xmin[i] = min( BH_x )
        BH_Xmax[i] = max( BH_x )
        BH_Ymin[i] = min( BH_y )
        BH_Ymax[i] = max( BH_y )
        if i==0:
            plt.plot( BH_x, BH_y, color='red',   label="BH"+str(i+1), linewidth=2 )
        elif i==1:
            plt.plot( BH_x, BH_y, color='green', label="BH"+str(i+1), linewidth=2 )
        elif i==2:
            plt.plot( BH_x, BH_y, color='blue',  label="BH"+str(i+1), linewidth=2 )
        elif i==3:
            plt.plot( BH_x, BH_y, color='gray',  label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "X [M]",          fontsize=16 )
    plt.ylabel( "Y [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Xmin0 = min( BH_Xmin )
    Xmax0 = max( BH_Xmax )
    Ymin0 = min( BH_Ymin )
    Ymax0 = max( BH_Ymax )
    Xmin  = min( Xmin0-2.0, -5.0 )
    Xmax  = max( Xmax0+2.0, +5.0 )
    Ymin  = min( Ymin0-2.0, -5.0 )
    Ymax  = max( Ymax0+2.0, +5.0 )
    plt.xlim( Xmin, Xmax )          # x 轴范围从 Xmin 到 Xmax
    plt.ylim( Ymin, Ymax )          # y 轴范围从 Ymin 到 Ymax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    # plt.show(                                                      )
    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_XY.pdf") )
    plt.close(                                                       )
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（XZ图）
    
    plt.figure( figsize=(8,8)                         )   ## 这里 figsize 可以设定图形的大小
    plt.title( " Black Hole Trajectory ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    
    for i in range(input_data.puncture_number):
        BH_x       = data[:, 3*i+1]
        BH_y       = data[:, 3*i+2]
        BH_z       = data[:, 3*i+3]
        BH_Xmin[i] = min( BH_x )
        BH_Xmax[i] = max( BH_x )
        BH_Zmin[i] = min( BH_z )
        BH_Zmax[i] = max( BH_z )
        if i==0:
            plt.plot( BH_x, BH_z, color='red',   label="BH"+str(i+1), linewidth=2 )
        elif i==1:
            plt.plot( BH_x, BH_z, color='green', label="BH"+str(i+1), linewidth=2 )
        elif i==2:
            plt.plot( BH_x, BH_z, color='blue',  label="BH"+str(i+1), linewidth=2 )
        elif i==3:
            plt.plot( BH_x, BH_z, color='gray',  label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "X [M]",          fontsize=16 )
    plt.ylabel( "Z [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Xmin0 = min( BH_Xmin )
    Xmax0 = max( BH_Xmax )
    Zmin0 = min( BH_Zmin )
    Zmax0 = max( BH_Zmax )
    Xmin  = min( Xmin0-2.0, -5.0 )
    Xmax  = max( Xmax0+2.0, +5.0 )
    Zmin  = min( Zmin0-2.0, -5.0 )
    Zmax  = max( Zmax0+2.0, +5.0 )
    plt.xlim( Xmin, Xmax )         # x 轴范围从 Xmin 到 Xmax
    plt.ylim( Zmin, Zmax )         # y 轴范围从 Zmin 到 Zmax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    # plt.show(                                                      )
    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_XZ.pdf") )
    plt.close(                                                       )
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（YZ图）
    
    plt.figure( figsize=(8,8)                         )   ## 这里 figsize 可以设定图形的大小
    plt.title( " Black Hole Trajectory ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    
    for i in range(input_data.puncture_number):
        BH_x       = data[:, 3*i+1]
        BH_y       = data[:, 3*i+2]
        BH_z       = data[:, 3*i+3]
        BH_Ymin[i] = min( BH_y )
        BH_Ymax[i] = max( BH_y )
        BH_Zmin[i] = min( BH_z )
        BH_Zmax[i] = max( BH_z )
        if i==0:
            plt.plot( BH_y, BH_z, color='red',   label="BH"+str(i+1), linewidth=2 )
        elif i==1:
            plt.plot( BH_y, BH_z, color='green', label="BH"+str(i+1), linewidth=2 )
        elif i==2:
            plt.plot( BH_y, BH_z, color='blue',  label="BH"+str(i+1), linewidth=2 )
        elif i==3:
            plt.plot( BH_y, BH_z, color='gray',  label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "Y [M]",          fontsize=16 )
    plt.ylabel( "Z [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Ymin0 = min( BH_Ymin )
    Ymax0 = max( BH_Ymax )
    Zmin0 = min( BH_Zmin )
    Zmax0 = max( BH_Zmax )
    Ymin  = min( Ymin0-2.0, -5.0 )
    Ymax  = max( Ymax0+2.0, +5.0 )
    Zmin  = min( Zmin0-2.0, -5.0 )
    Zmax  = max( Zmax0+2.0, +5.0 )
    plt.xlim( Ymin, Ymax )          # x 轴范围从 Ymin 到 Ymax
    plt.ylim( Zmin, Zmax )          # y 轴范围从 Zmin 到 Zmax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    # plt.show(                                                      )
    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_YZ.pdf") )
    plt.close(                                                       )
    
    # --------------------------
    
    # 得到黑洞 1 和黑洞 2 的坐标
    BH_x1 = data[:, 1]
    BH_y1 = data[:, 2]
    BH_z1 = data[:, 3]
    BH_x2 = data[:, 4]
    BH_y2 = data[:, 5]
    BH_z2 = data[:, 6]
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（X2-X1 Y2-Y1）

    plt.figure( figsize=(8,8)                                           )                          
    plt.title(  " Black Hole Trajectory ",                  fontsize=18 )   
    plt.plot(   (BH_x2-BH_x1), (BH_y2-BH_y1), color='blue', linewidth=2 )
    plt.xlabel( " $X_{2}$ - $X_{1}$ [M] ",                  fontsize=16 )
    plt.ylabel( " $Y_{2}$ - $Y_{1}$ [M] ",                  fontsize=16 )
    plt.legend( loc='upper right'                                       )

    # 设置坐标轴的范围
    Xmin0 = min( (BH_x2 - BH_x1) )
    Xmax0 = max( (BH_x2 - BH_x1) ) 
    Ymin0 = min( (BH_y2 - BH_y1) )
    Ymax0 = max( (BH_y2 - BH_y1) ) 
    Xmin  = min( Xmin0-2.0, -5.0 )
    Xmax  = max( Xmax0+2.0, +5.0 )
    Ymin  = min( Ymin0-2.0, -5.0 )
    Ymax  = max( Ymax0+2.0, +5.0 )
    plt.xlim( Xmin, Xmax )          # x 轴范围从 Xmin 到 Xmax
    plt.ylim( Ymin, Ymax )          # y 轴范围从 Ymin 到 Ymax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_21_XY.pdf")  )
    plt.close(                                                           )
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（X2-X1 Z2-Z1）
    
    plt.figure( figsize=(8,8)                                           )                          
    plt.title(  " Black Hole Trajectory ",                  fontsize=18 )   
    plt.plot(   (BH_x2-BH_x1), (BH_z2-BH_z1), color='blue', linewidth=2 )
    plt.xlabel( " $X_{2}$ - $X_{1}$ [M] ",                  fontsize=16 )
    plt.ylabel( " $Z_{2}$ - $Z_{1}$ [M] ",                  fontsize=16 )
    plt.legend( loc='upper right'                                       )

    # 设置坐标轴的范围
    Xmin0 = min( (BH_x2 - BH_x1) )
    Xmax0 = max( (BH_x2 - BH_x1) ) 
    Zmin0 = min( (BH_z2 - BH_z1) )
    Zmax0 = max( (BH_z2 - BH_z1) ) 
    Xmin  = min( Xmin0-2.0, -5.0 )
    Xmax  = max( Xmax0+2.0, +5.0 )
    Zmin  = min( Zmin0-2.0, -5.0 )
    Zmax  = max( Zmax0+2.0, +5.0 )
    plt.xlim( Xmin, Xmax )          # x 轴范围从 Xmin 到 Xmax
    plt.ylim( Zmin, Zmax )          # y 轴范围从 Zmin 到 Zmax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_21_XZ.pdf")  )
    plt.close(                                                           )
    
    # --------------------------
    
    # 画出黑洞位移的轨迹图（Y2-Y1 Z2-Z1）
    
    plt.figure( figsize=(8,8)                                           )                          
    plt.title(  " Black Hole Trajectory ",                  fontsize=18 )   
    plt.plot(   (BH_y2-BH_y1), (BH_z2-BH_z1), color='blue', linewidth=2 )
    plt.xlabel( " $Y_{2}$ - $Y_{1}$ [M] ",                  fontsize=16 )
    plt.ylabel( " $Z_{2}$ - $Z_{1}$ [M] ",                  fontsize=16 )
    plt.legend( loc='upper right'                                       )

    # 设置坐标轴的范围
    Ymin0 = min( (BH_y2 - BH_y1) )
    Ymax0 = max( (BH_y2 - BH_y1) ) 
    Zmin0 = min( (BH_z2 - BH_z1) )
    Zmax0 = max( (BH_z2 - BH_z1) ) 
    Ymin  = min( Ymin0-2.0, -5.0 )
    Ymax  = max( Ymax0+2.0, +5.0 )
    Zmin  = min( Zmin0-2.0, -5.0 )
    Zmax  = max( Zmax0+2.0, +5.0 )
    plt.xlim( Ymin, Ymax )          # x 轴范围从 Ymin 到 Ymax
    plt.ylim( Zmin, Zmax )          # y 轴范围从 Zmin 到 Zmax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_21_YZ.pdf")  )
    plt.close(                                                           )
    
    # --------------------------
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    print(                      )
    print( " 对黑洞轨迹画图完成 " )
    print( " Black holes' trajectory plot has been finished (2D plot)" )
    print(                                                             )

    return

####################################################################################



####################################################################################

## 该函数对黑洞的相对距离画图

def generate_puncture_distence_plot( outdir, figure_outdir ):

    print(                                                )
    print( " 正在对黑洞间距进行画图 "                        )
    print( " Plotting the black hole relative distance " )
    print(                                               )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    print( " 对应数据文件为 ",              file0 )
    print( " Corresponding data file = ", file0 )

    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)
    
    # --------------------------
    
    # 画出每个黑洞距坐标原点的距离 R 随时间的变化图

    # 初始化黑洞距离的最大值和最小值
    BH_Rmin = numpy.zeros(input_data.puncture_number)
    BH_Rmax = numpy.zeros(input_data.puncture_number)
    
    # 创建一个新的图
    fig = plt.figure( figsize=(8,8) )
    plt.title( " Black Hole Position R ", fontsize=18 )   # 添加标题
    
    BH_time = data[:, 0]
    
    for i in range(input_data.puncture_number):
        BH_x = data[:, 3*i+1]
        BH_y = data[:, 3*i+2]
        BH_z = data[:, 3*i+3]
        BH_R = (BH_x*BH_x + BH_y*BH_y + BH_z*BH_z)**0.5
        # 利用 numpy 直接平方和求出距离 R
        BH_Rmin[i] = min( BH_R )
        BH_Rmax[i] = max( BH_R )
        if i==0:
            plt.plot( BH_time, BH_R, color='red',   label="BH"+str(i+1), linewidth=2 )
        elif i==1:
            plt.plot( BH_time, BH_R, color='green', label="BH"+str(i+1), linewidth=2 )
        elif i==2:
            plt.plot( BH_time, BH_R, color='blue',  label="BH"+str(i+1), linewidth=2 )
        elif i==3:
            plt.plot( BH_time, BH_R, color='gray',  label="BH"+str(i+1), linewidth=2 )

    # 设置坐标轴标签
    plt.xlabel( " $T$ [M] ",      fontsize=16 )
    plt.ylabel( " $R$ [M] ",      fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    R_min0 = min( BH_Rmin ) 
    R_max0 = max( BH_Rmax )
    R_min  = max( R_min0-2.0,  0.0 )
    R_max  = max( R_max0+2.0, +5.0 )
    plt.ylim( R_min, R_max )             # y 轴范围从 R_min 到 R_max
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    # plt.show(                                                   )
    plt.savefig( os.path.join(figure_outdir, "BH_Position_R.pdf") )
    plt.close(                                                    )
    
    # --------------------------
    
    # 得到黑洞 1 和黑洞 2 的坐标
    BH_x1  = data[:, 1]
    BH_y1  = data[:, 2]
    BH_z1  = data[:, 3]
    BH_x2  = data[:, 4]
    BH_y2  = data[:, 5]
    BH_z2  = data[:, 6]
    
    # 利用 numpy 直接计算平方和开根号，得到黑洞 1 和黑洞 2 的相对距离
    BH_R12 = ( (BH_x2-BH_x1)**2 + (BH_y2-BH_y1)**2 + (BH_z2*BH_z1)**2 )**0.5
    
    # --------------------------
    
    # 画出黑洞 1 和黑洞 2 的相对轨迹图 R12

    plt.figure( figsize=(8,8)                              )                          
    plt.title(  " Black Hole Distance ",       fontsize=18 )   
    plt.plot(   BH_time, BH_R12, color='blue', linewidth=2 )
    plt.xlabel( " $T$ [M] ",                   fontsize=16 )
    plt.ylabel( " $R_{12}$ [M] ",              fontsize=16 )
    plt.legend( loc='upper right'                          )

    # 设置坐标轴的范围
    R12_min0 = min( BH_R12 )
    R12_max0 = max( BH_R12 ) 
    R12_min  = max( R12_min0-2.0,  0.0 )
    R12_max  = max( R12_max0+2.0, +5.0 )
    plt.ylim( R12_min, R12_max )             # y 轴范围从 R12_min 到 R12_max
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Distance_21.pdf")  )
    plt.close(                                                      )
    
    print(                           )
    print( " 正在对黑洞间距画图完成 "   )
    print( " black hole relative distance plot has been finished " )
    print(                                                         )
    
    # --------------------------
 
    return

####################################################################################



####################################################################################

## 该函数对黑洞轨迹画 3 维图

def generate_puncture_orbit_plot3D( outdir, figure_outdir ):

    print(                               )
    print( " 正在对黑洞轨迹进行画 3 维图 " )
    print( " Plotting the black holes' trajectory (3D plot) " )
    print(                               )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    print(" 对应数据文件为 ",               file0)
    print( " Corresponding data file = ", file0 )

    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)

    # 初始化黑洞坐标的最大值和最小值
    BH_Xmin = numpy.zeros(input_data.puncture_number)
    BH_Xmax = numpy.zeros(input_data.puncture_number)
    BH_Ymin = numpy.zeros(input_data.puncture_number)
    BH_Ymax = numpy.zeros(input_data.puncture_number)
    BH_Zmin = numpy.zeros(input_data.puncture_number)
    BH_Zmax = numpy.zeros(input_data.puncture_number)
    
    # 创建一个新的图
    fig = plt.figure( figsize=(8,8) )
 
    # 创建一个3D坐标轴
    ax = fig.add_subplot(111, projection='3d')
    # 添加标题
    ax.set_title( " Black Hole Trajectory ", fontsize=18 )
    
    for i in range(input_data.puncture_number):
        BH_x = data[:, 3*i+1]
        BH_y = data[:, 3*i+2]
        BH_z = data[:, 3*i+3]
        BH_Xmin[i] = min( BH_x )
        BH_Xmax[i] = max( BH_x )
        BH_Ymin[i] = min( BH_y )
        BH_Ymax[i] = max( BH_y )
        BH_Zmin[i] = min( BH_z )
        BH_Zmax[i] = max( BH_z )
        if i==0:
            ax.plot( BH_x, BH_y, BH_z, color='red',   label="BH"+str(i+1), linewidth=2 )
        elif i==1:
            ax.plot( BH_x, BH_y, BH_z, color='green', label="BH"+str(i+1), linewidth=2 )
        elif i==2:
            ax.plot( BH_x, BH_y, BH_z, color='blue',  label="BH"+str(i+1), linewidth=2 )
        elif i==3:
            ax.plot( BH_x, BH_y, BH_z, color='gray',  label="BH"+str(i+1), linewidth=2 )

    # 设置坐标轴标签
    ax.set_xlabel( "X [M]",          fontsize=16 )
    ax.set_ylabel( "Y [M]",          fontsize=16 )
    ax.set_zlabel( "Z [M]",          fontsize=16 )
    plt.legend(    loc='upper right'             )

    # 设置坐标轴的范围
    Xmin0 = min( BH_Xmin )
    Xmax0 = max( BH_Xmax )
    Ymin0 = min( BH_Ymin )
    Ymax0 = max( BH_Ymax )
    Zmin0 = min( BH_Zmin )
    Zmax0 = max( BH_Zmax )
    Xmin  = min( Xmin0-2.0, -5.0 )
    Xmax  = max( Xmax0+2.0, +5.0 )
    Ymin  = min( Ymin0-2.0, -5.0 )
    Ymax  = max( Ymax0+2.0, +5.0 )
    Zmin  = min( Zmin0-2.0, -5.0 )
    Zmax  = max( Zmax0+2.0, +5.0 )
    ax.set_xlim( [Xmin, Xmax] )      # y 轴范围从 Ymin 到 Ymax
    ax.set_ylim( [Ymin, Ymax] )      # y 轴范围从 Ymin 到 Ymax
    ax.set_zlim( [Zmin, Zmax] )      # z 轴范围从 Zmin 到 Zmax

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_3D.pdf") )
    plt.close(                                                       )
    
    print(                                                             )
    print( " 对黑洞轨迹 3 维画图完成 "                                    )
    print( " Black holes' trajectory plot has been finished (3D plot)" )
    print(                                                             )
 
    return


####################################################################################



####################################################################################

## 该函数对引力波波形 Psi4 画图

def generate_gravitational_wave_psi4_plot( outdir, figure_outdir, detector_number_i ):
    

    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_psi4.dat")

    if ( detector_number_i == 0 ):
        print(                                                )
        print( " 对 Weyl 共形变量 Psi4 进行画图 "               )
        print( " Plotting the Weyl conformal component Psi4 " )
        print(                                                )
        print( " 对应数据文件为 ",             file0 )
        print( " corresponding data file = ", file0 )
        print(                                      )

    print( " 对第 ", detector_number_i, " 个探测器半径数据进行画图 " )
    print( " Begin the Weyl conformal Psi4 plot for detector number = ", detector_number_i )
    
    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)
    
    # 取出 phi4 文件中各列的数据
    time                 = data[:,0]
    psi4_l2m2m_real      = data[:,1]
    psi4_l2m2m_imaginary = data[:,2]
    psi4_l2m1m_real      = data[:,3]
    psi4_l2m1m_imaginary = data[:,4]
    psi4_l2m0_real       = data[:,5]
    psi4_l2m0_imaginary  = data[:,6]
    psi4_l2m1_real       = data[:,7]
    psi4_l2m1_imaginary  = data[:,8]
    psi4_l2m2_real       = data[:,9]
    psi4_l2m2_imaginary  = data[:,10]
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    # python 中除法会返回浮点数，因此这里设置为整除
    length = len(time) // input_data.Detector_Number 
    
    time2                 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2m_real2      = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2m_imaginary2 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1m_real2      = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1m_imaginary2 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m0_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m0_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    
    # 将数据拆分为各探测器半径对应的数据
    for i in range(input_data.Detector_Number):
        for j in range(length):
            time2[i,j]                 = time[                 j*input_data.Detector_Number + i ]
            psi4_l2m2m_real2[i,j]      = psi4_l2m2m_real[      j*input_data.Detector_Number + i ]
            psi4_l2m2m_imaginary2[i,j] = psi4_l2m2m_imaginary[ j*input_data.Detector_Number + i ]
            psi4_l2m1m_real2[i,j]      = psi4_l2m1m_real[      j*input_data.Detector_Number + i ]
            psi4_l2m1m_imaginary2[i,j] = psi4_l2m1m_imaginary[ j*input_data.Detector_Number + i ]
            psi4_l2m0_real2[i,j]       = psi4_l2m0_real[       j*input_data.Detector_Number + i ]
            psi4_l2m0_imaginary2[i,j]  = psi4_l2m0_imaginary[  j*input_data.Detector_Number + i ]
            psi4_l2m1_real2[i,j]       = psi4_l2m1_real[       j*input_data.Detector_Number + i ]
            psi4_l2m1_imaginary2[i,j]  = psi4_l2m1_imaginary[  j*input_data.Detector_Number + i ]
            psi4_l2m2_real2[i,j]       = psi4_l2m2_real[       j*input_data.Detector_Number + i ]
            psi4_l2m2_imaginary2[i,j]  = psi4_l2m2_imaginary[  j*input_data.Detector_Number + i ]
            
    # 根据输入数据推算出探测器距离
    Detector_Interval   = ( input_data.Detector_Rmax - input_data.Detector_Rmin ) / ( input_data.Detector_Number - 1 )
    Detector_Distance_R = input_data.Detector_Rmax - Detector_Interval * detector_number_i
    
    plt.figure( figsize=(8,8) )                                   ## 这里 figsize 可以设定图形的大小
    plt.title( f" Gravitational Wave $\Psi_{4}$   Detector Distance =  { Detector_Distance_R } ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    plt.plot( time2[detector_number_i], psi4_l2m0_real2[detector_number_i],      \
              color='red',    label="l=2 m=0 real",                       linewidth=2 )
    plt.plot( time2[detector_number_i], psi4_l2m0_imaginary2[detector_number_i], \
              color='orange', label="l=2 m=0 imaginary",  linestyle='--', linewidth=2 )
    plt.plot( time2[detector_number_i], psi4_l2m1_real2[detector_number_i],      \
              color='green',  label="l=2 m=1 real",                       linewidth=2 )
    plt.plot( time2[detector_number_i], psi4_l2m1_imaginary2[detector_number_i], \
              color='cyan',   label="l=2 m=1 imaginary",  linestyle='--', linewidth=2 )
    plt.plot( time2[detector_number_i], psi4_l2m2_real2[detector_number_i],      \
              color='black',  label="l=2 m=2 real",                       linewidth=2 )
    plt.plot( time2[detector_number_i], psi4_l2m2_imaginary2[detector_number_i], \
              color='gray',   label="l=2 m=2 imaginary",  linestyle='--', linewidth=2 )
    plt.xlabel( "T [M]",          fontsize=16 )
    plt.ylabel( r"$R*\Psi$",      fontsize=16 )
    plt.legend( loc='upper right'             )
    plt.grid(   color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线
    plt.savefig( os.path.join(figure_outdir, "Gravitational_Psi4_Detector_" + str(detector_number_i) + ".pdf") )
    
    
    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " The Weyl Conformal component Psi4 plot has been finished ", " detector number ", detector_number_i )
    print(                                                                                            )

    if ( detector_number_i == (input_data.Detector_Number-1) ):
        print(                                    )
        print( " 对 Weyl 共形变量 Psi4 的画图完成 " )
        print( " The Weyl conformal component Psi4 plots have been finished " )
        print(                                                                 )


    return

####################################################################################



####################################################################################

## 该函数对时空 ADM 质量画图

def generate_ADMmass_plot( outdir, figure_outdir, detector_number_i ):

    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_ADMQs.dat")

    if ( detector_number_i == 0 ):
        print(                                                )
        print( " 对时空 ADM 质量和角动量进行画图 "              )
        print( " Plotting the ADM mass and angular momentum " )
        print(                                                )
        print( " 对应数据文件为 ",             file0 )
        print( " corresponding data file = ", file0 )
        print(                                      )
    
    print( " 对第 ", detector_number_i, " 个探测器半径数据进行画图 " )
    print( " Begin the ADM momentum plot for detector number =  ", detector_number_i )


    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)
    
    # 取出 AMD 动量文件中各列的数据
    time     = data[:,0]
    ADM_mass = data[:,1]
    ADM_Px   = data[:,2]
    ADM_Py   = data[:,3]
    ADM_Pz   = data[:,4]
    ADM_Jx   = data[:,5]
    ADM_Jy   = data[:,6]
    ADM_Jz   = data[:,7]
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    # python 中除法会返回浮点数，因此这里设置为整除
    length = len(time) // input_data.Detector_Number
    
    '''
    # 将数据拆分为各探测器半径对应的数据
    time2     = time.reshape( (input_data.Detector_Number, length) )
    ADM_mass2 = ADM_mass.reshape( (input_data.Detector_Number, length) )
    ADM_Px2   = ADM_Px.reshape( (input_data.Detector_Number, length) )
    ADM_Py2   = ADM_Py.reshape( (input_data.Detector_Number, length) )
    ADM_Pz2   = ADM_Pz.reshape( (input_data.Detector_Number, length) )
    ADM_Jx2   = ADM_Jx.reshape( (input_data.Detector_Number, length) )
    ADM_Jy2   = ADM_Jy.reshape( (input_data.Detector_Number, length) )
    ADM_Jz2   = ADM_Jz.reshape( (input_data.Detector_Number, length) )
    '''
    # reshape 的行和列没有搞清楚，换成笨办法
    time2     = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_mass2 = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Px2   = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Py2   = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Pz2   = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Jx2   = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Jy2   = numpy.zeros( (input_data.Detector_Number, length) )
    ADM_Jz2   = numpy.zeros( (input_data.Detector_Number, length) )
    
    # 将数据拆分为各探测器半径对应的数据
    for i in range(input_data.Detector_Number):
        for j in range(length):
            time2[i,j]     = time[     j*input_data.Detector_Number + i ]
            ADM_mass2[i,j] = ADM_mass[ j*input_data.Detector_Number + i ]
            ADM_Px2[i,j]   = ADM_Px[   j*input_data.Detector_Number + i ]
            ADM_Py2[i,j]   = ADM_Py[   j*input_data.Detector_Number + i ]
            ADM_Pz2[i,j]   = ADM_Pz[   j*input_data.Detector_Number + i ]
            ADM_Jx2[i,j]   = ADM_Jx[   j*input_data.Detector_Number + i ]
            ADM_Jy2[i,j]   = ADM_Jy[   j*input_data.Detector_Number + i ]
            ADM_Jz2[i,j]   = ADM_Jz[   j*input_data.Detector_Number + i ]
            
    # 根据输入数据推算出探测器距离
    Detector_Interval   = ( input_data.Detector_Rmax - input_data.Detector_Rmin ) / ( input_data.Detector_Number - 1 )
    Detector_Distance_R = input_data.Detector_Rmax - Detector_Interval * detector_number_i
            
    # 画出当前 detector 半径的 AMD 动量
    plt.figure( figsize=(8,8) )                  
    plt.title(f" ADM Momentum    Detector Distence = {Detector_Distance_R}", fontsize=18 )   
    plt.plot( time2[detector_number_i], ADM_mass2[detector_number_i], color='red',   label="ADM Mass", linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Px2[detector_number_i],   color='green', label="ADM Px",   linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Py2[detector_number_i],   color='cyan',  label="ADM Py",   linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Pz2[detector_number_i],   color='blue',  label="ADM Pz",   linewidth=2 )
    plt.xlabel( "T [M]",            fontsize=16 )
    plt.ylabel( "ADM Momentum [M]", fontsize=16 )
    plt.legend( loc='upper right'               )
    plt.grid(   color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线
    plt.savefig( os.path.join(figure_outdir, "ADM_Mass_Dector_" + str(detector_number_i) + ".pdf") )
    
    # 画出当前 detector 半径的 AMD 角动量
    plt.figure( figsize=(8,8) )                  
    plt.title(f" ADM Angular Momentum    Detector Distence = {Detector_Distance_R}", fontsize=18 )   
    # plt.plot( time2[detector_number_i], ADM_mass2[detector_number_i], color='red',   label="ADM Mass", linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Jx2[detector_number_i],   color='green', label="ADM Jx",   linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Jy2[detector_number_i],   color='cyan',  label="ADM Jy",   linewidth=2 )
    plt.plot( time2[detector_number_i], ADM_Jz2[detector_number_i],   color='blue',  label="ADM Jz",   linewidth=2 )
    plt.xlabel( "T [M]",                        fontsize=16 )
    plt.ylabel( "ADM Angular Momentum [$M^2$]", fontsize=16 )
    plt.legend( loc='upper right'                           )
    plt.grid(   color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线
    plt.savefig( os.path.join(figure_outdir, "ADM_Angular_Momentum_Dector_" + str(detector_number_i) + ".pdf") )
    

    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " ADM momentum plot has been finished, detector number =  ", detector_number_i )
    print(                                                                                )

    if ( detector_number_i == (input_data.Detector_Number-1) ):
        print( " 对时空 ADM 质量和角动量的画图完成 " )
        print( " The ADM mass and augular momentum plots have been finished " )
        print(                                                                )

    return
    
####################################################################################



####################################################################################

## 该函数对哈密顿约束违反性况画图

def generate_constraint_check_plot( outdir, figure_outdir, input_level_number ):

    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_constraint.dat")

    if ( input_level_number == 0 ):
        print(                                                   )
        print( " 对哈密顿约束违反性况进行画图 "                     )
        print( " Plotting the constraint violation for each grid level" )
        print(                                                          )
        print( " 对应数据文件为 ",             file0 )
        print( " corresponding data file = ", file0 )
        print(                                      )

    print( " 对第 ", input_level_number, " 层网格数据进行画图 " )
    print( " Begin the constraint violation plot for grid level number =  ", input_level_number )
    
    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)
    
    # 取出约束数据文件中各列的数据
    time          = data[:,0]
    Constraint_H  = data[:,1]
    Constraint_Px = data[:,2]
    Constraint_Py = data[:,3]
    Constraint_Pz = data[:,4]
    Constraint_Gx = data[:,5]
    Constraint_Gy = data[:,6]
    Constraint_Gz = data[:,7]
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    # 初始化各类数据
    
    if (input_data.basic_grid_set == "Patch"):
        level_number = input_level_number   
        length0      = input_data.grid_level
        # python 中除法会返回浮点数，因此这里设置为整除
        length1      = len(time) // length0   
    elif (input_data.basic_grid_set == "Shell-Patch"):
        # 如果格点类型选择为 Shell-Patch，网格层数加 1
        level_number = input_level_number + 1
        length0      = input_data.grid_level + 1 
        # python 中除法会返回浮点数，因此这里设置为整除
        length1      = len(time) // length0   
    
    time2          = numpy.zeros( (length0, length1) )
    Constraint_H2  = numpy.zeros( (length0, length1) )
    Constraint_Px2 = numpy.zeros( (length0, length1) )
    Constraint_Py2 = numpy.zeros( (length0, length1) )
    Constraint_Pz2 = numpy.zeros( (length0, length1) )
    Constraint_Gx2 = numpy.zeros( (length0, length1) )
    Constraint_Gy2 = numpy.zeros( (length0, length1) )
    Constraint_Gz2 = numpy.zeros( (length0, length1) )
    
    # 将数据拆分为各探测器半径对应的数据
    for i in range(length0):
        for j in range(length1):
            time2[i,j]          = time[          j*length0 + i ]
            Constraint_H2[i,j]  = Constraint_H[  j*length0 + i ]
            Constraint_Px2[i,j] = Constraint_Px[ j*length0 + i ]
            Constraint_Py2[i,j] = Constraint_Py[ j*length0 + i ]
            Constraint_Pz2[i,j] = Constraint_Pz[ j*length0 + i ]
    
    # 画出最外 detector 半径的约束违反
    plt.figure( figsize=(8,8) )                    
    plt.title( f" ADM Constraint  Grid Level = {input_level_number}", fontsize=18 )   
    plt.plot( time2[level_number], Constraint_H2[level_number],  color='red',   label="ADM Constraint H",  linewidth=2 )
    plt.plot( time2[level_number], Constraint_Px2[level_number], color='green', label="ADM Constraint Px", linewidth=2 )
    plt.plot( time2[level_number], Constraint_Py2[level_number], color='cyan',  label="ADM Constraint Py", linewidth=2 )
    plt.plot( time2[level_number], Constraint_Pz2[level_number], color='blue',  label="ADM Constraint Pz", linewidth=2 )
    plt.xlabel( "T [M]",          fontsize=16 )
    plt.ylabel( "ADM Constraint", fontsize=16 )
    plt.legend( loc='upper right'             )
    plt.grid(   color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线
    plt.savefig( os.path.join(figure_outdir, "ADM_Constraint_Grid_Level_" + str(input_level_number) + ".pdf") )
    

    print( " 第 ", input_level_number, " 层网格数据的约束违反情况的画图完成 " )
    print( " Constraint violation plot has been finished, grid level number = ", input_level_number )
    print(                                                                                          )
    
    if ( input_level_number == (input_data.grid_level-1) ):
        print( " 对约束违反情况的画图完成 " )
        print( " Constraint violation plot has been finished " )
        print(                                                 )

    return

####################################################################################



####################################################################################

# 单独使用的例子
'''
outdir = "./BBH_q=1"

generate_puncture_orbit_plot(    outdir, outdir )
generate_puncture_orbit_plot3D(  outdir, outdir )
generate_puncture_distence_plot( outdir, outdir )

for i in range(input_data.grid_level):
    generate_constraint_check_plot( outdir, outdir, i )

for i in range(input_data.Detector_Number):
    generate_ADMmass_plot( outdir, outdir, i )

for i in range(input_data.Detector_Number):
    generate_gravitational_wave_psi4_plot( outdir, outdir, i )
'''
####################################################################################


