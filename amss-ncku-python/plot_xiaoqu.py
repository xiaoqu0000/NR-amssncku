
#################################################
##
## 这个文件对 AMSS-NCKU 数值相对论的结果进行画图
## 小曲
## 2024/10/01 --- 2026/02/27
## 2026/08/21 修改为2进制并行画图
##
#################################################

import numpy                               ## 导入 numpy 包进行数组的操作
import matplotlib                          ## 导入 matplotlib 包
matplotlib.use('Agg')                      ## 强制使用无界面的 Agg 后端，避免 GUI 后端（如 gtk4agg）在多进程 fork 画图时死锁
import matplotlib.pyplot    as     plt     ## 导入 matplotlib 包进行画图
from   mpl_toolkits.mplot3d import Axes3D  ## 画 3 维图需要
import glob
import os                                  ## 导入 os 包进行系统操作
import multiprocessing                     ## 导入 multiprocessing 包进行多进程并行画图

import plot_binary_data
import AMSS_NCKU_Input as input_data

# plt.rcParams['text.usetex'] = True  ## 在绘图中允许使用 latex 字体



####################################################################################

## 该函数是在并行画图中被子进程调用的函数
## 每个子进程负责一个二进制文件的画图任务

def _plot_binary_data_worker( args ):

    input_language, filename, binary_outdir, figure_outdir = args

    plot_binary_data.plot_binary_data( input_language, filename, binary_outdir, figure_outdir )

    return filename

####################################################################################

## 该函数对若干二进制文件的画图任务进行多进程并行处理
## 并行画图不改变每个文件的画图算法和输出内容

def _plot_binary_data_parallel( tasks ):

    ## 没有任务时直接返回
    if ( len(tasks) == 0 ):
        return

    ## 进程数目取输入文件中设定的进程数与任务数目的较小值
    if ( input_data.plot_binary_data_processes > 0 ):
        nproc = min( input_data.plot_binary_data_processes, len(tasks) )
    else:
        ## 如果输入文件中进程数目设为 0，则自动使用所有 CPU 核数
        nproc = min( os.cpu_count(), len(tasks) )

    ## 任务数目不超过 1 个时，直接按顺序画图
    if ( nproc <= 1 ):
        for task in tasks:
            _plot_binary_data_worker( task )
        return

    if ( tasks[0][0] == "Chinese" ):
        print(                                             )
        print( " 使用 ", nproc, " 个进程并行画图 "          )
        print(                                             )
    else:
        print(                                                  )
        print( " plotting in parallel with ", nproc, " processes " )
        print(                                                  )

    ## 利用进程池并行处理所有画图任务
    with multiprocessing.Pool( processes=nproc ) as pool:
        for _ in pool.imap_unordered( _plot_binary_data_worker, tasks, chunksize=1 ):
            pass

    return

####################################################################################

## 该函数根据二进制数据画出所有二维图

def generate_binary_data_plot( input_language, binary_outdir, figure_outdir ):


    if ( input_language == "Chinese" ):
        print(                                        )
        print( " 开始 AMSS-NCKU 程序输出的二进制数据画图 " )
        print(                                        )
    elif ( input_language == "English" ):
        print(                                                              )
        print( " Beginning the AMSS-NCKU Binary Data Plotting From Output " )
        print(                                                              )

    ###########################################

    ## 设定画图的文件夹目录
    
    ## 利用 numpy 设置空数组，用于存放字符串
    ## 设置 dtype='U200'，确保存放的是字符串，长度为 200（防止长度不够）
    figure_outdir_static  = numpy.empty( input_data.static_grid_level, dtype='U200' )
    figure_outdir_moving  = numpy.empty( input_data.moving_grid_level, dtype='U200' )
    figure_outdir_moving2 = numpy.empty( (input_data.moving_grid_level, input_data.puncture_number), dtype='U200' )
    
    ## 设置每层网格画图文件夹的目录
    for i in range(input_data.static_grid_level):
        figure_outdir_static[i] = os.path.join( figure_outdir, "level" + str(i) )
    for i in range(input_data.moving_grid_level):
        j = i + input_data.static_grid_level
        for k in range(input_data.puncture_number):
            figure_outdir_moving[i] = os.path.join( figure_outdir, "level" + str(j) )
            figure_outdir_moving2[i,k] = os.path.join( figure_outdir_moving[i], "puncture" + str(k) )

    ###########################################

    if ( input_language == "Chinese" ):
        print(                                  )
        print( " 读取 AMSS-NCKU 程序的二进制数据 " )
        print(                                  )
    elif ( input_language == "English" ):
        print(                                               )
        print( " Reading AMSS-NCKU Binary Data From Output " )
        print(                                               )
    
    if ( input_language == "Chinese" ):
        print(                               ) 
        print( " 输出文件中的二进制数据文件列表 " )
        print(                               )
    elif ( input_language == "English" ):
        print(                                        )
        print( " The output binary data files list: " )
        print(                                        )

    ## 初始化文件列表
    file_list_all   = []
    ## 生成用于存放字符串的空数组
    file_list_static_grid = [ [] for _ in range(input_data.static_grid_level) ]
    ## empty_list = [[] for _ in range(m)]  # 创建5个空列表的列表
    file_list_moving_grid = [] 
    for i in range(input_data.moving_grid_level):
        xx = [ [] for _ in range(input_data.puncture_number) ]
        file_list_moving_grid.append(xx)
    
    ## 抓取所有的文件名称
    globby_all_files = glob.glob( os.path.join(binary_outdir, '*.bin') ) 

    ## 将抓取到的文件名输出，并组成列表
    file_list_all = []
    for filename in sorted(globby_all_files):
        filename1 = filename.replace(binary_outdir + "/", "") # 去掉路径中的前缀
        file_list_all.append(filename1)
        print(filename1)

    ## 测试代码
    ## print( globby_all_files )
    ## print( file_list_all    )

    for i in range(input_data.static_grid_level):

        if ( input_language == "Chinese" ):
            print(                                           ) 
            print( " 输出文件中的二进制数据文件列表  对应网格层", i )
            print(                                           ) 
        elif ( input_language == "English" ):
            print(                                                       ) 
            print( " The output binary data files list in grid level", i )
            print(                                                       ) 

        ## 筛选出该层固定网格对应的二进制数据文件
        for filename in sorted(file_list_all):
            ## 筛选开头字母为 “Lev-0x” 的文件
            ## 利用 fstring 来补充确保整数的长度为 2
            if filename.startswith( "Lev" + f"{i:02d}" ):
                file_list_static_grid[i].append(filename)    
                print(filename)

    for i in range(input_data.moving_grid_level):
                
        j = i + input_data.static_grid_level

        for k in range(input_data.puncture_number):

            if ( input_language == "Chinese" ):
                print(                                                           )
                print( " 输出文件中的二进制数据文件列表  对应网格层", j, " puncture", k )
                print(                                                           )
            elif ( input_language == "English" ):
                print(                                                                       )
                print( " The output binary data files list in grid level", j, " puncture", k )
                print(                                                                       )

            ## 筛选出该层移动网格对应的二进制数据文件
            for filename in sorted(file_list_all):
                if filename.startswith( "Lev" + f"{j:02d}" + "-" + f"{k:02d}" ):
                    file_list_moving_grid[i][k].append(filename)   
                    print(filename)
    
    ## 使用 glob 抓取文件的例子
    '''
    # 第一步：抓取当前目录下所有 .bin 文件
    bin_files = glob.glob("*.bin")

    # 第二步：在第一步结果基础上，筛选文件名以 "abc" 开头的文件
    # 使用 str.startswith() 方法精确匹配前缀
    abc_bin_files = [f for f in bin_files if f.lower().startswith("abc")]
    '''

    ###########################################

    ## 收集所有需要画图的二进制文件任务
    tasks = []

    ## 如果输入文件中 plot_binary_data_level 设定为 "All-Level"，对所有网格层的二进制文件进行画图

    if (input_data.plot_binary_data_level == "All-Level"): 
        
        ## 对固定网格层的二进制数画图
        for i in range(input_data.static_grid_level):

            ## 生成该层二进制数据画图的文件夹目录
            os.mkdir( figure_outdir_static[i] )
            ## 生成该层二进制数据画图的各种子文件夹目录
            plot_binary_data.generate_binary_data_plot_directionary( figure_outdir_static[i] )

            ## 收集该网格层的所有二进制文件的画图任务
            for filename in file_list_static_grid[i]:     
                tasks.append( ( input_language, filename, binary_outdir, figure_outdir_static[i] ) )

        ## 对移动网格层的二进制数画图
        for i in range(input_data.moving_grid_level):

            ## 生成该层二进制数据画图的文件夹目录
            os.mkdir( figure_outdir_moving[i] )
            
            j = i + input_data.static_grid_level
            
            for k in range(input_data.puncture_number):

                ## 生成该层二进制数据画图的各种子文件夹目录
                os.mkdir( figure_outdir_moving2[i,k] )
                plot_binary_data.generate_binary_data_plot_directionary( figure_outdir_moving2[i,k] )

                ## 收集该网格层的所有二进制文件的画图任务
                for filename in file_list_moving_grid[i][k]:     
                    tasks.append( ( input_language, filename, binary_outdir, figure_outdir_moving2[i,k] ) )

    ###########################################

    ## 如果输入文件中 plot_binary_data_level 设定为 "Single-Level"，对设定的网格层的二进制文件进行画图

    elif (input_data.plot_binary_data_level == "Single-Level"):

        ## 设定的网格层在固定网格层
        if (input_data.plot_binary_data_levelnumber < input_data.static_grid_level):

            i = input_data.plot_binary_data_levelnumber

            ## 生成该层二进制数据画图的文件夹目录
            os.mkdir( figure_outdir_static[i] )
            ## 生成该层二进制数据画图的各种子文件夹目录
            plot_binary_data.generate_binary_data_plot_directionary( figure_outdir_static[i] )

            ## 收集该网格层的所有二进制文件的画图任务
            for filename in file_list_static_grid[i]:     
                tasks.append( ( input_language, filename, binary_outdir, figure_outdir_static[i] ) )

        ## 设定的网格层在移动网格层
        else:
            
            j = input_data.plot_binary_data_levelnumber - input_data.static_grid_level

            ## 生成该层二进制数据画图的文件夹目录
            os.mkdir( figure_outdir_moving[j] )
            
            for k in range(input_data.puncture_number):

                ## 生成该层二进制数据画图的各种子文件夹目录
                os.mkdir( figure_outdir_moving2[j,k] )
                plot_binary_data.generate_binary_data_plot_directionary( figure_outdir_moving2[j,k] )

                ## 收集该网格层的所有二进制文件的画图任务
                for filename in file_list_moving_grid[j][k]:     
                    tasks.append( ( input_language, filename, binary_outdir, figure_outdir_moving2[j,k] ) )

    ###########################################

    ## 对所有收集到的二进制文件画图任务进行多进程并行画图
    _plot_binary_data_parallel( tasks )

    if ( input_language == "Chinese" ):
        print(                         )
        print( " 二进制数据画图已完成 " )
        print(                         )
    elif ( input_language == "English" ):
        print(                                        )
        print( " Binary Data Plot Has been Finished " )
        print(                                        )

    return

####################################################################################



####################################################################################

## 该函数对黑洞轨迹画图

def generate_puncture_orbit_plot( input_language, outdir, figure_outdir ):

    if ( input_language == "Chinese" ):
        print(                        )
        print( " 正在对黑洞轨迹进行画图 " )
        print(                        )
    elif ( input_language == "English" ):
        print(                                                   )
        print( " Plotting the black holes' trajectory (2D plot)" )
        print(                                                   )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    if ( input_language == "Chinese" ):
        print( " 对应数据文件为 ",              file0 )
    elif ( input_language == "English" ):
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

    # 设定曲线的颜色
    line_color = [ 'red', 'green', 'blue', 'black', 'gray', 'yellow', 'cyan', 'pink', 'magenta' ]
    
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
        plt.plot( BH_x, BH_y, color=line_color[i], label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "X [M]",          fontsize=16 )
    plt.ylabel( "Y [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Xmin0 = min( BH_Xmin )
    Xmax0 = max( BH_Xmax )
    Ymin0 = min( BH_Ymin )
    Ymax0 = max( BH_Ymax )
    Xmin  = min( Xmin0-2.0, -3.0 )
    Xmax  = max( Xmax0+2.0, +3.0 )
    Ymin  = min( Ymin0-2.0, -3.0 )
    Ymax  = max( Ymax0+2.0, +3.0 )
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
        plt.plot( BH_x, BH_z, color=line_color[i], label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "X [M]",          fontsize=16 )
    plt.ylabel( "Z [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Xmin0 = min( BH_Xmin )
    Xmax0 = max( BH_Xmax )
    Zmin0 = min( BH_Zmin )
    Zmax0 = max( BH_Zmax )
    Xmin  = min( Xmin0-2.0, -3.0 )
    Xmax  = max( Xmax0+2.0, +3.0 )
    Zmin  = min( Zmin0-2.0, -3.0 )
    Zmax  = max( Zmax0+2.0, +3.0 )
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
        plt.plot( BH_y, BH_z, color=line_color[i], label="BH"+str(i+1), linewidth=2 )
            
    plt.xlabel( "Y [M]",          fontsize=16 )
    plt.ylabel( "Z [M]",          fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    Ymin0 = min( BH_Ymin )
    Ymax0 = max( BH_Ymax )
    Zmin0 = min( BH_Zmin )
    Zmax0 = max( BH_Zmax )
    Ymin  = min( Ymin0-2.0, -3.0 )
    Ymax  = max( Ymax0+2.0, +3.0 )
    Zmin  = min( Zmin0-2.0, -3.0 )
    Zmax  = max( Zmax0+2.0, +3.0 )
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
    Xmin  = min( Xmin0-2.0, -3.0 )
    Xmax  = max( Xmax0+2.0, +3.0 )
    Ymin  = min( Ymin0-2.0, -3.0 )
    Ymax  = max( Ymax0+2.0, +3.0 )
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
    Xmin  = min( Xmin0-2.0, -3.0 )
    Xmax  = max( Xmax0+2.0, +3.0 )
    Zmin  = min( Zmin0-2.0, -3.0 )
    Zmax  = max( Zmax0+2.0, +3.0 )
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
    Ymin  = min( Ymin0-2.0, -3.0 )
    Ymax  = max( Ymax0+2.0, +3.0 )
    Zmin  = min( Zmin0-2.0, -3.0 )
    Zmax  = max( Zmax0+2.0, +3.0 )
    plt.xlim( Ymin, Ymax )          # x 轴范围从 Ymin 到 Ymax
    plt.ylim( Zmin, Zmax )          # y 轴范围从 Zmin 到 Zmax
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_21_YZ.pdf")  )
    plt.close(                                                           )
    
    # --------------------------
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    if ( input_language == "Chinese" ):
        print(                       )
        print( " 对黑洞轨迹画图完成 " )
        print(                       )
    elif ( input_language == "English" ):
        print(                                                             )
        print( " Black holes' trajectory plot has been finished (2D plot)" )
        print(                                                             )

    return

####################################################################################



####################################################################################

## 该函数对黑洞的相对距离画图

def generate_puncture_distence_plot( input_language, outdir, figure_outdir ):

    if ( input_language == "Chinese" ):
        print(                           )
        print( " 正在对黑洞间距进行画图 " )
        print(                           )
    elif ( input_language == "English" ):
        print(                                               )
        print( " Plotting the black hole relative distance " )
        print(                                               )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    if ( input_language == "Chinese" ):
        print( " 对应数据文件为 ",             file0 )
    elif ( input_language == "English" ):
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

    # 设定曲线的颜色
    line_color = [ 'red', 'green', 'blue', 'black', 'gray', 'yellow', 'cyan', 'pink', 'magenta' ]
    
    for i in range(input_data.puncture_number):
        BH_x = data[:, 3*i+1]
        BH_y = data[:, 3*i+2]
        BH_z = data[:, 3*i+3]
        BH_R = (BH_x*BH_x + BH_y*BH_y + BH_z*BH_z)**0.5
        # 利用 numpy 直接平方和求出距离 R
        BH_Rmin[i] = min( BH_R )
        BH_Rmax[i] = max( BH_R )
        plt.plot( BH_time, BH_R, color=line_color[i], label="BH"+str(i+1), linewidth=2 )
        

    # 设置坐标轴标签
    plt.xlabel( " $T$ [M] ",      fontsize=16 )
    plt.ylabel( " $R$ [M] ",      fontsize=16 )
    plt.legend( loc='upper right'             )

    # 设置坐标轴的范围
    R_min0 = min( BH_Rmin ) 
    R_max0 = max( BH_Rmax )
    R_min  = max( R_min0-2.0,  0.0 )
    R_max  = max( R_max0+2.0, +3.0 )
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
    R12_max  = max( R12_max0+2.0, +3.0 )
    plt.ylim( R12_min, R12_max )             # y 轴范围从 R12_min 到 R12_max
    
    plt.grid( color='gray', linestyle='--', linewidth=0.5 )  # 显示网格线

    plt.savefig( os.path.join(figure_outdir, "BH_Distance_21.pdf")  )
    plt.close(                                                      )
    
    if ( input_language == "Chinese" ):
        print(                           )
        print( " 正在对黑洞间距画图完成 " )
        print(                           )
    elif ( input_language == "English" ):
        print(                                                         )
        print( " black hole relative distance plot has been finished " )
        print(                                                         )
    
    # --------------------------
 
    return

####################################################################################



####################################################################################

## 该函数对黑洞轨迹画 3 维图

def generate_puncture_orbit_plot3D( input_language, outdir, figure_outdir ):

    if ( input_language == "Chinese" ):
        print(                             )
        print( " 正在对黑洞轨迹进行画 3 维图 " )
        print(                             )
    elif ( input_language == "English" ):
        print(                                                    )
        print( " Plotting the black holes' trajectory (3D plot) " )
        print(                                                    )
    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_BH.dat")
    
    if ( input_language == "Chinese" ):
        print(" 对应数据文件为 ",               file0 )
    elif ( input_language == "English" ):
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
    ax = fig.add_subplot( 111, projection='3d' )
    # 添加标题
    ax.set_title( " Black Hole Trajectory ", fontsize=18 )

    # 设定曲线的颜色
    line_color = [ 'red', 'green', 'blue', 'black', 'gray', 'yellow', 'cyan', 'pink', 'magenta' ]
    
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
        ax.plot( BH_x, BH_y, BH_z, color=line_color[i], label="BH"+str(i+1), linewidth=2 )
        

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
    Xmin  = min( Xmin0-2.0, -3.0 )
    Xmax  = max( Xmax0+2.0, +3.0 )
    Ymin  = min( Ymin0-2.0, -3.0 )
    Ymax  = max( Ymax0+2.0, +3.0 )
    Zmin  = min( Zmin0-2.0, -3.0 )
    Zmax  = max( Zmax0+2.0, +3.0 )
    ax.set_xlim( [Xmin, Xmax] )      # y 轴范围从 Ymin 到 Ymax
    ax.set_ylim( [Ymin, Ymax] )      # y 轴范围从 Ymin 到 Ymax
    ax.set_zlim( [Zmin, Zmax] )      # z 轴范围从 Zmin 到 Zmax

    plt.savefig( os.path.join(figure_outdir, "BH_Trajectory_3D.pdf") )
    plt.close(                                                       )
    
    if ( input_language == "Chinese" ):
        print(                          )
        print( " 对黑洞轨迹 3 维画图完成 " )
        print(                          )
    elif ( input_language == "English" ):
        print(                                                             )
        print( " Black holes' trajectory plot has been finished (3D plot)" )
        print(                                                             )
 
    return


####################################################################################



####################################################################################

## 该函数对引力波波形 Psi4 画图

def generate_gravitational_wave_psi4_plot( input_language, outdir, figure_outdir, detector_number_i ):
    

    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_psi4.dat")

    if ( detector_number_i == 0 ):
        if ( input_language == "Chinese" ):
            print(                                  )
            print( " 对 Weyl 共形变量 Psi4 进行画图 " )
            print(                                  )
            print( " 对应数据文件为 ", file0         )
            print(                                  )
        elif ( input_language == "English" ):
            print(                                                )
            print( " Plotting the Weyl conformal component Psi4 " )
            print(                                                )
            print( " corresponding data file = ", file0           )
            print(                                                )

    if ( input_language == "Chinese" ):
        print( " 对第 ", detector_number_i, " 个探测器半径数据进行画图 " )
    elif ( input_language == "English" ):
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
    
    
    if ( input_language == "Chinese" ):
        print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
        print(                                                     )
    elif ( input_language == "English" ):
        print( " The Weyl Conformal component Psi4 plot has been finished ", " detector number ", detector_number_i )
        print(                                                                                                      )

    if ( detector_number_i == (input_data.Detector_Number-1) ):
        if ( input_language == "Chinese" ):
            print(                                   )
            print( " 对 Weyl 共形变量 Psi4 的画图完成 " )
            print(                                   )
        elif ( input_language == "English" ):
            print(                                                                )
            print( " The Weyl conformal component Psi4 plots have been finished " )
            print(                                                                )


    return

####################################################################################



####################################################################################

## 该函数对时空 ADM 质量画图

def generate_ADMmass_plot( input_language, outdir, figure_outdir, detector_number_i ):

    
    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_ADMQs.dat")

    if ( detector_number_i == 0 ):
        if ( input_language == "Chinese" ):
            print(                                 )
            print( " 对时空 ADM 质量和角动量进行画图 " )
            print(                                 )
            print( " 对应数据文件为 ", file0         )
            print(                                 )
        elif ( input_language == "English" ):
            print(                                                )
            print( " Plotting the ADM mass and angular momentum " )
            print(                                                )
            print( " corresponding data file = ", file0           )
            print(                                                )
    
    if ( input_language == "Chinese" ):
        print( " 对第 ", detector_number_i, " 个探测器半径数据进行画图 " )
    elif ( input_language == "English" ):
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
    

    if ( input_language == "Chinese" ):
        print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
        print(                                                        )
    elif ( input_language == "English" ):
        print( " ADM momentum plot has been finished, detector number =  ", detector_number_i )
        print(                                                                                )

    if ( detector_number_i == (input_data.Detector_Number-1) ):
        if ( input_language == "Chinese" ):
            print( " 对时空 ADM 质量和角动量的画图完成 " )
            print(                                  )
        elif ( input_language == "English" ):
            print( " The ADM mass and augular momentum plots have been finished " )
            print(                                                                )

    return
    
####################################################################################



####################################################################################

## 该函数对哈密顿约束违反性况画图

def generate_constraint_check_plot( input_language, outdir, figure_outdir, input_level_number ):

    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_constraint.dat")

    if ( input_level_number == 0 ):
        if ( input_language == "Chinese" ):
            print(                             )
            print( " 对哈密顿约束违反性况进行画图 " )
            print(                             )
            print( " 对应数据文件为 ", file0     )
            print(                             )
        elif ( input_language == "English" ):
            print(                                                          )
            print( " Plotting the constraint violation for each grid level" )
            print(                                                          )
            print( " corresponding data file = ", file0                     )
            print(                                                          )

    if ( input_language == "Chinese" ):
        print( " 对第 ", input_level_number, " 层网格数据进行画图 " )
    elif ( input_language == "English" ):
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
    

    if ( input_language == "Chinese" ):
        print( " 第 ", input_level_number, " 层网格数据的约束违反情况的画图完成 " )
        print(                                                                  )
    elif ( input_language == "English" ):
        print( " Constraint violation plot has been finished, grid level number = ", input_level_number )
        print(                                                                                          )
    
    if ( input_level_number == (input_data.grid_level-1) ):
        if ( input_language == "Chinese" ):
            print( " 对约束违反情况的画图完成 " )
            print(                          )
        elif ( input_language == "English" ):
            print( " Constraint violation plot has been finished " )
            print(                                                 )

    return

####################################################################################



####################################################################################

# 单独使用的例子
'''
outdir = "./GW190521-Hierarchical/AMSS_NCKU_output/"
figure_dir = "./GW190521-Hierarchical/figure/"

input_language = 'English'

generate_puncture_orbit_plot(    input_language, outdir, figure_dir )
generate_puncture_orbit_plot3D(  input_language, outdir, figure_dir )
generate_puncture_distence_plot( input_language, outdir, figure_dir )

for i in range(input_data.grid_level):
    generate_constraint_check_plot( input_language, outdir, figure_dir, i )

for i in range(input_data.Detector_Number):
    generate_ADMmass_plot( input_language, outdir, figure_dir, i )

for i in range(input_data.Detector_Number):
    generate_gravitational_wave_psi4_plot( input_language, outdir, figure_dir, i )
'''
####################################################################################


