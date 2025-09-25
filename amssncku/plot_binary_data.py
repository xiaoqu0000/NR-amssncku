
#################################################
##
## 这个文件包含了数值相对论所的二进制数据画图
## 小曲
## 2024/10/01 --- 2025/09/14 
##
#################################################

import numpy
import scipy
import matplotlib.pyplot    as     plt
from   matplotlib.colors    import LogNorm
from   mpl_toolkits.mplot3d import Axes3D
## import torch
import AMSS_NCKU_Input      as input_data

import os


#########################################################################################

def plot_binary_data( filename, binary_outdir, figure_outdir ):

    figure_title0 = filename.replace(binary_outdir + "/", "")  # 去掉路径中的前缀
    figure_title  = figure_title0.replace(".bin", "")          # 去掉路径中的.bin
    
    print(                                                    )
    print( " 正在读取二进制文件 = ",               figure_title0 )
    print( " reading binary data from file = ", figure_title0 )

###################################

    # 打开文件
    # 根据 AMSS-NCKU 输出二进制文件中的数据顺序依次读入数据
    with open(filename, 'rb') as file:

        physical_time = numpy.fromfile( file, dtype=numpy.float64, count=1 )
        nx, ny, nz    = numpy.fromfile( file, dtype=numpy.int32,   count=3 )
        xmin, xmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        ymin, ymax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        zmin, zmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        data          = numpy.fromfile( file, dtype=numpy.float64          )
        
        # 现在 data 数组包含了文件中的二进制数据
 
    print( " 读取的数组大小 = ",     data.shape                            ) 
    print( " 读取的数组长度 = ",     data.size                             ) 
    print( " 原始设定的数组长度 = ", nx, "*", ny, "*", nz, " = ", nx*ny*nz )
    print( " obtained data shape  = ",     data.shape                        ) 
    print( " obtained data size   = ",     data.size                         ) 
    print( " obtained data points = ", nx, "*", ny, "*", nz, " = ", nx*ny*nz )
    
###################################

    # 将读入的数据转化为多维数组
    data_reshape = data.reshape( (nz, ny, nx) ) ## 经过测试，这样的排列方式画出来才正常（reshape的第一个需要是z方向）
    # print(data_reshape)

    # data1 = data_reshape[0,:,:]
    # print(data1)

    Rmin = [xmin, ymin, zmin] 
    Rmax = [xmax, ymax, zmax]
    N    = [nx, ny, nz]
    print( " 格点坐标最小值 = ", Rmin  )
    print( " 格点坐标最大值 = ", Rmax  )
    print( " 格点数目       = ", N    )
    print( " coordinate minimum = ", Rmin )
    print( " coordinate maximum = ", Rmax )
    print( " grid point         = ", N    )
    
    print(                                                )
    print( " 数据读取完成，接下来开始画图 "                   )
    print( " Data file read successfully. Plotting data " )
    print(                                                )

    # 利用画图函数进行画图
    figure_title0    = filename.replace(binary_outdir + "/", "") # 去掉路径中的前缀
    figure_title     = figure_title.replace(".bin", "")          # 去掉最后的".bin"
    figure_title_new = figure_title[:-6]                         # 再去掉末尾的6个字符，代表的是迭代次数
    
    get_data_xy( Rmin, Rmax, N, data_reshape, physical_time[0], figure_title_new, figure_outdir )
    # 注意 numpy 从二进制文件中读取的 physical_time 是一个数组（尽管实际上只有一个元素）
    # 因此用 physical_time[0] 代表对应的时间值
    
    # 手动删除数据以清除内存
    del data
    del data_reshape
    
    print( " 二进制文件 = ", figure_title0," 画图已完成 "                 )
    print( " binary data file = ", figure_title0," plot has finished " )
    print(                                                             )

    return
    
    
#########################################################################################




####################################################################################

# 这是一个对某一二进制数据的画图函数

def get_data_xy( Rmin, Rmax, n, data0, time, figure_title, figure_outdir ):

    figure_contourplot_outdir = os.path.join(figure_outdir, "contour plot")
    figure_densityplot_outdir = os.path.join(figure_outdir, "density plot")
    figure_surfaceplot_outdir = os.path.join(figure_outdir, "surface plot")

    # 根据读到的格点信息还原格点坐标
    x = numpy.linspace(Rmin[0], Rmax[0], n[0])
    y = numpy.linspace(Rmin[1], Rmax[1], n[1])
    z = numpy.linspace(Rmin[2], Rmax[2], n[2])
    # print(x)
    # print(y)
    # print(z)

    # 用 meshgrid 建立二维格点坐标                             
    # X, Y = torch.meshgrid(torch.tensor(x), torch.tensor(y))    # 除了 numpy 以外，torch 也可以 meshgrid
    X, Y = numpy.meshgrid(x, y)    
    
    # 补充 numpy.meshgrid 的相关信息
    # 假设 x 是长度为 nx 的数组，y 为长度为 ny 的数组
    # 而 X, Y = numpy.meshgrid(x, y) 得到的 X 和 Y 都是 (ny, nx) 数组，并且 X 的每一行都是 x 的副本，Y 的每一列都是 y 的副本   
    
    # print( X.shape )
    # print( Y.shape )
    # print( X[0,:] )
    # print( Y[:,0] )

    # 获取 xy 平面上的数据
    if input_data.Symmetry == "no-symmetry":
        data_xy = data0[n[2]//2,:,:]
    else:
        data_xy = data0[0,:,:]
        
    # 由于原始数据是按照列排列的，因此我们转置之后再画图
    # 经过测试后发现无需转置
    # data_xy = numpy.transpose(data_xy_0)
    
    # print( data_xy_0.shape )
    # print( data_xy.shape )
    
    # 设置新坐标，方便对数据进行插值 
    x_new = numpy.linspace(Rmin[0], Rmax[0], int(2.5*n[0]))
    y_new = numpy.linspace(Rmin[1], Rmax[1], int(2.5*n[1]))
    z_new = numpy.linspace(Rmin[2], Rmax[2], int(2.5*n[2]))
    X_new, Y_new = numpy.meshgrid(x_new, y_new)
    
    # 对数据进行插值
    data_xy_fit = scipy.interpolate.griddata( (X.flatten(), Y.flatten()), data_xy.flatten(), (X_new, Y_new), method="cubic" )

    # 下面画出二维等高线图
    fig, ax = plt.subplots()
    # contourf = ax.contourf(X, Y, data_xy, 8, cmap='coolwarm', norm=LogNorm(vmin=1, vmax=10), levels=numpy.logspace(-2, 2, 8))  # 使用'coolwarm'色板，并设置标准色彩映射
    # contourf = ax.contourf( X, Y, data_xy_0, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Y, data_xy_0, 8, colors='k', linewidths=0.5 )     # 添加等高线
    # 由于原始数据是按照列排列的，因此我们转置之后再画图
    # contourf = ax.contourf( X, Y, data_xy, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Y, data_xy, 8, colors='k', linewidths=0.5 )       # 添加等高线
    # 对插值后的数据进行画图
    contourf = ax.contourf( X_new, Y_new, data_xy_fit, cmap=plt.get_cmap('RdYlGn_r') )
    contour  = ax.contour(  X_new, Y_new, data_xy_fit, 8, colors='k', linewidths=0.5 )     # 添加等高线
    cbar     = plt.colorbar(contourf)                                                      # 添加色条
    ax.set_title(  figure_title + "  physical time = " + str(time) )                       # 设置标题和轴标签
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    # plt.show()                                                                                                          # 显示图像
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " contour_plot.pdf") )   # 保存图像
    plt.close()
    
    # 下面画出二维热图    
    # fig1 = plt.figure()
    fig1, ax  = plt.subplots()
    # 经过测试后发现，似乎不用转置，直接用 imshowfig 画图就可以得到相应的 xy 图
    # imshowfig = plt.imshow( data_xy, interpolation='bicubic', extent=[X.min(), X.max(), Y.min(), Y.max()] )
    # imshowfig = plt.imshow( numpy.transpose(data_xy), interpolation='bicubic', extent=[X.min(), X.max(), Y.min(), Y.max()] )
    # 不知道为什么，y轴坐标是反的，需要人为处理以下才能化成正确的图
    imshowfig = plt.imshow( data_xy, interpolation='bicubic', extent=[X.min(), X.max(), Y.max(), Y.min()] )                                                    
    ax.invert_yaxis()                                                       # 将 y 轴的正负翻转
    cbar      = plt.colorbar(imshowfig)                                     # 添加色条
    ax.set_title(  figure_title + "  physical time = " + str(time)  )       # 设置标题和轴标签
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    # plt.show() 
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " density_plot.pdf") )
    plt.close()

    # 下面画出三维图
    fig2 = plt.figure()                                                       # 创建一个新的图像
    ax = fig2.add_subplot( 111, projection='3d' )                             # 创建一个 3D 绘图区域
    # 对插值后的数据进行画图
    # ax.plot_surface( X, Y, data_xy, cmap='viridis' )                        # 绘制曲面
    ax.plot_surface( X_new, Y_new, data_xy_fit, cmap='viridis' )              # 绘制曲面
    ax.set_title(  figure_title + "  physical time = " + str(time) )          # 设置标题和轴标签
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    # plt.show()                                                              # 显示图像
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " surface_plot.pdf") )   # 保存图像
    plt.close()

    return

####################################################################################

