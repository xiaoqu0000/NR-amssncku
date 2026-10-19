
#################################################
##
## 这个文件包含了数值相对论所的二进制数据画图 | This file contains the binary-data plotting needed for numerical relativity
## 小曲 | Xiao Qu
## 2024/10/01 --- 2026/03/10 
##
#################################################

import numpy
import scipy
import matplotlib                                       ## 导入 matplotlib 包 | Import the matplotlib package
matplotlib.use('Agg')                                   ## 强制使用无界面的 Agg 后端，避免 GUI 后端（如 gtk4agg）在多进程 fork 画图时死锁 | Force the headless Agg backend to avoid deadlocks of GUI backends (e.g. gtk4agg) when plotting with multiprocessing fork
import matplotlib.pyplot    as     plt
from   matplotlib.colors    import LogNorm
from   mpl_toolkits.mplot3d import Axes3D
## import torch
import AMSS_NCKU_Input      as input_data

import os


#########################################################################################

def plot_binary_data( input_language, filename, binary_outdir, figure_outdir ):

    
    if ( input_language == "Chinese" ):
        print(                                 )
        print( " 正在读取二进制文件 = ", filename )
    elif ( input_language == "English" ):
        print(                                               )
        print( " reading binary data from file = ", filename )


    ## 由于函数传入的文件名不包含路径中前缀，需要补充前缀来读取文件 | Since the filename passed into the function does not contain the path prefix, the prefix must be added back to read the file
    filename_oringinal = os.path.join( binary_outdir, filename )

    ###################################

    # 打开文件 | Open the file
    # 根据 AMSS-NCKU 输出二进制文件中的数据顺序依次读入数据 | Read in the data one by one according to the data order in the AMSS-NCKU binary output file
    with open(filename_oringinal, 'rb') as file:

        physical_time = numpy.fromfile( file, dtype=numpy.float64, count=1 )
        nx, ny, nz    = numpy.fromfile( file, dtype=numpy.int32,   count=3 )
        xmin, xmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        ymin, ymax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        zmin, zmax    = numpy.fromfile( file, dtype=numpy.float64, count=2 )
        data          = numpy.fromfile( file, dtype=numpy.float64          )
        
        # 现在 data 数组包含了文件中的二进制数据 | Now the data array contains the binary data from the file

        # === 转为 Python 原生数值 === | Convert to native Python numeric values
        physical_time = physical_time.item()
        xmin, xmax = xmin.item(), xmax.item()
        ymin, ymax = ymin.item(), ymax.item()
        zmin, zmax = zmin.item(), zmax.item()
        nx, ny, nz = int(nx), int(ny), int(nz)

 
    if ( input_language == "Chinese" ):
        print( " 读取的数组大小 = ",    data.shape                            ) 
        print( " 读取的数组长度 = ",    data.size                             ) 
        print( " 原始设定的数组长度 = ", nx, "*", ny, "*", nz, " = ", nx*ny*nz )
    elif ( input_language == "English" ):
        print( " obtained data shape  = ", data.shape                            ) 
        print( " obtained data size   = ", data.size                             ) 
        print( " obtained data points = ", nx, "*", ny, "*", nz, " = ", nx*ny*nz )
    
    ###################################

    # 将读入的数据转化为多维数组 | Reshape the read-in data into a multidimensional array
    data_reshape = data.reshape( (nz, ny, nx) ) ## 经过测试，这样的排列方式画出来才正常（reshape的第一个需要是z方向） | Tests show the figure looks right only with this ordering (the first argument of reshape must be the z direction)
    # print(data_reshape)

    # data1 = data_reshape[0,:,:]
    # print(data1)

    Rmin = [xmin, ymin, zmin] 
    Rmax = [xmax, ymax, zmax]
    N    = [nx, ny, nz]
    
    if ( input_language == "Chinese" ):
        print( " 格点坐标最小值 = ", Rmin  )
        print( " 格点坐标最大值 = ", Rmax  )
        print( " 格点数目       = ", N    )
    elif ( input_language == "English" ):
        print( " coordinate minimum = ", Rmin )
        print( " coordinate maximum = ", Rmax )
        print( " grid point         = ", N    )
    
    if ( input_language == "Chinese" ):
        print(                             )
        print( " 数据读取完成，接下来开始画图 " )
        print(                             )
    elif ( input_language == "English" ):
        print(                                                )
        print( " Data file read successfully. Plotting data " )
        print(                                                )

    ###################################

    ## 利用画图函数进行画图 | Use the plotting functions to draw the figures

    ## 设定文件名称 | Set the file name
    ## 经过修改后，该函数传入的文件名已经不包含路径中前缀 | After the modification, the filename passed into this function no longer contains the path prefix
    ## figure_title0 = filename.replace(binary_outdir + "/", "")  # 去掉路径中的前缀 | Remove the prefix from the path
    ## figure_title  = figure_title0.replace(".bin", "")          # 去掉最后的".bin" | Remove the trailing ".bin"
    
    figure_title     = filename.replace(".bin", "")              # 去掉路径中的.bin | Remove ".bin" from the path
    figure_title1    = figure_title[:-6]                         # 再去掉末尾的6个字符，代表的是迭代次数 | Then remove the last 6 characters, which represent the iteration number
    figure_title_new = figure_title1[9:]                         # 再去掉开头的9个字符，代表的是"Lev0x-0x_" | Then remove the first 9 characters, which represent "Lev0x-0x_"

    ## 设定文件夹目录 | Set the folder directories
    figure_outdir_xy = os.path.join( figure_outdir, "XY plot" )
    figure_outdir_xz = os.path.join( figure_outdir, "XZ plot" )
    figure_outdir_yz = os.path.join( figure_outdir, "YZ plot" )
    
    ## 根据输入文件设定对哪个平面的数据进行画图 | According to the input file, decide which plane of the data to plot
    if (input_data.plot_binary_data_set == "xy-xz-yz-plot"):
        get_data_xy( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_xy )
        get_data_xz( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_xz )
        get_data_yz( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_yz )
    elif (input_data.plot_binary_data_set == "xy-plot"):
        get_data_xy( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_xy )
    elif (input_data.plot_binary_data_set == "xz-plot"):
        get_data_xz( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_xz )
    elif (input_data.plot_binary_data_set == "yz-plot"):
        get_data_yz( Rmin, Rmax, N, data_reshape, physical_time, figure_title_new, figure_outdir_yz )
    ## 注意前面已经将读取的 physical_time 转换为实数 | Note that the physical_time read in above has already been converted to a real number
    ## 因此不再需要 physical_time[0] 来代表对应的时间值 | So physical_time[0] is no longer needed to represent the corresponding time value

    ###################################
    
    # 手动删除数据以清除内存 | Manually delete the data to free the memory
    del data
    del data_reshape
    
    if ( input_language == "Chinese" ):
        print( " 二进制文件 = ", filename," 画图已完成 " )
        print(                                        )
    elif ( input_language == "English" ):
        print( " binary data file = ", filename," plot has finished " )
        print(                                                        )

    return
    
    
#########################################################################################




####################################################################################

## 这是一个对某一二进制数据的画图函数 | This is a plotting function for a given binary data file
## 对 xy 平面数据进行画图 | Plot the data on the xy plane

def get_data_xy( Rmin, Rmax, n, data0, time, figure_title, figure_outdir ):

    figure_contourplot_outdir = os.path.join(figure_outdir, "contour plot")
    figure_densityplot_outdir = os.path.join(figure_outdir, "density plot")
    figure_surfaceplot_outdir = os.path.join(figure_outdir, "surface plot")

    # 根据读到的格点信息还原格点坐标 | Reconstruct the grid coordinates from the read-in grid information
    x = numpy.linspace(Rmin[0], Rmax[0], n[0])
    y = numpy.linspace(Rmin[1], Rmax[1], n[1])
    z = numpy.linspace(Rmin[2], Rmax[2], n[2])
    # print(x)
    # print(y)
    # print(z)

    # 用 meshgrid 建立二维格点坐标 | Build the 2D grid coordinates with meshgrid                             
    # X, Y = torch.meshgrid(torch.tensor(x), torch.tensor(y))    # 除了 numpy 以外，torch 也可以 meshgrid | Besides numpy, torch can also do meshgrid
    X, Y = numpy.meshgrid(x, y)    
    
    # 补充 numpy.meshgrid 的相关信息 | Supplementary information on numpy.meshgrid
    # 假设 x 是长度为 nx 的数组，y 为长度为 ny 的数组 | Suppose x is an array of length nx and y is an array of length ny
    # 而 X, Y = numpy.meshgrid(x, y) 得到的 X 和 Y 都是 (ny, nx) 数组，并且 X 的每一行都是 x 的副本，Y 的每一列都是 y 的副本 | Then X, Y = numpy.meshgrid(x, y) gives X and Y both as (ny, nx) arrays, where each row of X is a copy of x and each column of Y is a copy of y   
    
    # print( X.shape )
    # print( Y.shape )
    # print( X[0,:] )
    # print( Y[:,0] )

    # 获取 xy 平面上的数据 | Get the data on the xy plane
    '''
    if input_data.Symmetry == "no-symmetry":
        data_xy = data0[n[2]//2,:,:]
    else:
        data_xy = data0[0,:,:]
    '''
    index = 0
    ## 找出最接近 z=0 的数据 | Find the data closest to z=0
    for i in range( n[2] ):
        if i > 0:
            if z[i-1] <= 0.0 < z[i]:
                if ( abs(z[i] - 0.0) <= abs(z[i-1] - 0.0) ):
                    index = i
                else:
                    index = i - 1
    data_xy = data0[index,:,:]
        
    ## 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    ## 经过测试后发现无需转置 | Tests showed that transposition is not needed
    ## data_xy_new = numpy.transpose(data_xy)
    
    ## print( data_xy.shape )
    ## print( data_xy_new.shape )
    
    ## 设置新坐标，方便对数据进行插值 | Set up new coordinates to make interpolating the data easier 
    x_new = numpy.linspace(Rmin[0], Rmax[0], int(2.5*n[0]))
    y_new = numpy.linspace(Rmin[1], Rmax[1], int(2.5*n[1]))
    z_new = numpy.linspace(Rmin[2], Rmax[2], int(2.5*n[2]))
    X_new, Y_new = numpy.meshgrid(x_new, y_new)
    
    ## 对数据进行插值 | Interpolate the data
    data_xy_fit = scipy.interpolate.griddata( (X.flatten(), Y.flatten()), data_xy.flatten(), (X_new, Y_new), method="cubic" )

    ## 下面画出二维等高线图 | Below, draw the 2D contour plot
    fig, ax = plt.subplots()
    # contourf = ax.contourf(X, Y, data_xy, 8, cmap='coolwarm', norm=LogNorm(vmin=1, vmax=10), levels=numpy.logspace(-2, 2, 8))  # 使用'coolwarm'色板，并设置标准色彩映射 | Use the 'coolwarm' colormap and set the standard color mapping
    # contourf = ax.contourf( X, Y, data_xy_0, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Y, data_xy_0, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    # 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    # contourf = ax.contourf( X, Y, data_xy, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Y, data_xy, 8, colors='k', linewidths=0.5 )       # 添加等高线 | Add contour lines
    # 对插值后的数据进行画图 | Plot the interpolated data
    contourf = ax.contourf( X_new, Y_new, data_xy_fit, cmap=plt.get_cmap('RdYlGn_r') )
    contour  = ax.contour(  X_new, Y_new, data_xy_fit, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    cbar     = plt.colorbar(contourf)                                                      # 添加色条 | Add a colorbar
    ax.set_title(  figure_title + "  physical time = " + str(time) )                       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " xy contour_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " xy contour_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()
    
    # 下面画出二维热图 | Below, draw the 2D heat map    
    # fig1 = plt.figure()
    fig1, ax  = plt.subplots()
    # 经过测试后发现，似乎不用转置，直接用 imshowfig 画图就可以得到相应的 xy 图 | Tests showed that transposition seems unnecessary: plotting directly with imshowfig gives the correct xy figure
    # imshowfig = plt.imshow( data_xy, interpolation='bicubic', extent=[X.min(), X.max(), Y.min(), Y.max()] )
    # imshowfig = plt.imshow( numpy.transpose(data_xy), interpolation='bicubic', extent=[X.min(), X.max(), Y.min(), Y.max()] )
    # 不知道为什么，y轴坐标是反的，需要人为处理以下才能化成正确的图 | For some reason the y-axis coordinates are reversed, and manual handling is needed to get a correct figure
    imshowfig = plt.imshow( data_xy, interpolation='bicubic', extent=[X.min(), X.max(), Y.max(), Y.min()] )                                                    
    cbar      = plt.colorbar(imshowfig)                                     # 添加色条 | Add a colorbar
    ax.invert_yaxis()                                                       # 将 y 轴的正负翻转 | Flip the sign of the y axis
    ax.set_title(  figure_title + "  physical time = " + str(time)  )       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " xy density_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " xy density_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    # 下面画出三维图 | Below, draw the 3D figure
    fig2 = plt.figure()                                                       # 创建一个新的图像 | Create a new figure
    ax = fig2.add_subplot( 111, projection='3d' )                             # 创建一个 3D 绘图区域 | Create a 3D plotting area
    # 对插值后的数据进行画图 | Plot the interpolated data
    # ax.plot_surface( X, Y, data_xy, cmap='viridis' )                        # 绘制曲面 | Draw the surface
    ax.plot_surface( X_new, Y_new, data_xy_fit, cmap='viridis' )              # 绘制曲面 | Draw the surface
    ax.set_title(  figure_title + "  physical time = " + str(time) )          # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Y [M]" )
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " xy surface_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " xy surface_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    return

####################################################################################



####################################################################################

## 这是一个对某一二进制数据的画图函数 | This is a plotting function for a given binary data file
## 对 xz 平面数据进行画图 | Plot the data on the xz plane

def get_data_xz( Rmin, Rmax, n, data0, time, figure_title, figure_outdir ):

    figure_contourplot_outdir = os.path.join(figure_outdir, "contour plot")
    figure_densityplot_outdir = os.path.join(figure_outdir, "density plot")
    figure_surfaceplot_outdir = os.path.join(figure_outdir, "surface plot")

    ## 根据读到的格点信息还原格点坐标 | Reconstruct the grid coordinates from the read-in grid information
    x = numpy.linspace(Rmin[0], Rmax[0], n[0])
    y = numpy.linspace(Rmin[1], Rmax[1], n[1])
    z = numpy.linspace(Rmin[2], Rmax[2], n[2])

    ## 用 meshgrid 建立二维格点坐标 | Build the 2D grid coordinates with meshgrid                             
    X, Z = numpy.meshgrid(x, z)    
    
    ## 补充 numpy.meshgrid 的相关信息 | Supplementary information on numpy.meshgrid
    ## 假设 x 是长度为 nx 的数组，y 为长度为 ny 的数组 | Suppose x is an array of length nx and y is an array of length ny
    ## 而 X, Y = numpy.meshgrid(x, y) 得到的 X 和 Y 都是 (ny, nx) 数组，并且 X 的每一行都是 x 的副本，Y 的每一列都是 y 的副本 | Then X, Y = numpy.meshgrid(x, y) gives X and Y both as (ny, nx) arrays, where each row of X is a copy of x and each column of Y is a copy of y   

    ## 获取 xz 平面上的数据 | Get the data on the xz plane
    '''
    data_xz = data0[:,0,:]
    '''
    index = 0
    ## 找出最接近 y=0 的数据 | Find the data closest to y=0
    for i in range( n[1] ):
        if i > 1:
            if y[i-1] <= 0.0 < y[i]:
                if ( abs(y[i] - 0.0) <= abs(y[i-1] - 0.0) ):
                    index = i
                else:
                    index = i - 1
    data_xz = data0[:,index,:]
                
        
    ## 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    ## 经过测试后发现无需转置 | Tests showed that transposition is not needed
    ## data_xz_new = numpy.transpose(data_xz)
    
    ## 设置新坐标，方便对数据进行插值 | Set up new coordinates to make interpolating the data easier 
    x_new = numpy.linspace(Rmin[0], Rmax[0], int(2.5*n[0]))
    y_new = numpy.linspace(Rmin[1], Rmax[1], int(2.5*n[1]))
    z_new = numpy.linspace(Rmin[2], Rmax[2], int(2.5*n[2]))
    X_new, Z_new = numpy.meshgrid(x_new, z_new)
    
    ## 对数据进行插值 | Interpolate the data
    data_xz_fit = scipy.interpolate.griddata( (X.flatten(), Z.flatten()), data_xz.flatten(), (X_new, Z_new), method="cubic" )

    ## 下面画出二维等高线图 | Below, draw the 2D contour plot
    fig, ax = plt.subplots()
    # contourf = ax.contourf( X, Z, data_xz, 8, cmap='coolwarm', norm=LogNorm(vmin=1, vmax=10), levels=numpy.logspace(-2, 2, 8))  # 使用'coolwarm'色板，并设置标准色彩映射 | Use the 'coolwarm' colormap and set the standard color mapping
    # contourf = ax.contourf( X, Z, data_xz_0, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Z, data_xz_0, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    # 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    # contourf = ax.contourf( X, Z, data_xz, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  X, Z, data_xz, 8, colors='k', linewidths=0.5 )       # 添加等高线 | Add contour lines
    # 对插值后的数据进行画图 | Plot the interpolated data
    contourf = ax.contourf( X_new, Z_new, data_xz_fit, cmap=plt.get_cmap('RdYlGn_r') )
    contour  = ax.contour(  X_new, Z_new, data_xz_fit, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    cbar     = plt.colorbar(contourf)                                                      # 添加色条 | Add a colorbar
    ax.set_title(  figure_title + "  physical time = " + str(time) )                       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " xz contour_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " xz contour_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()
    
    # 下面画出二维热图 | Below, draw the 2D heat map    
    # fig1 = plt.figure()
    fig1, ax  = plt.subplots()
    # 经过测试后发现，似乎不用转置，直接用 imshowfig 画图就可以得到相应的 xz 图 | Tests showed that transposition seems unnecessary: plotting directly with imshowfig gives the correct xz figure
    # imshowfig = plt.imshow( data_xz, interpolation='bicubic', extent=[X.min(), X.max(), Z.min(), Z.max()] )
    # imshowfig = plt.imshow( numpy.transpose(data_xz), interpolation='bicubic', extent=[X.min(), X.max(), Z.min(), Z.max()] )
    # 不知道为什么，z轴坐标是反的，需要人为处理以下才能化成正确的图 | For some reason the z-axis coordinates are reversed, and manual handling is needed to get a correct figure
    imshowfig = plt.imshow( data_xz, interpolation='bicubic', extent=[X.min(), X.max(), Z.max(), Z.min()] )                                                    
    cbar      = plt.colorbar(imshowfig)                                     # 添加色条 | Add a colorbar
    ax.invert_yaxis()                                                       # 将 y 轴的正负翻转（对应 z 的数据） | Flip the sign of the y axis (for the z data)
    ax.set_title(  figure_title + "  physical time = " + str(time)  )       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " xz density_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " xz density_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    # 下面画出三维图 | Below, draw the 3D figure
    fig2 = plt.figure()                                                       # 创建一个新的图像 | Create a new figure
    ax = fig2.add_subplot( 111, projection='3d' )                             # 创建一个 3D 绘图区域 | Create a 3D plotting area
    # 对插值后的数据进行画图 | Plot the interpolated data
    # ax.plot_surface( X, Z, data_xz, cmap='viridis' )                        # 绘制曲面 | Draw the surface
    ax.plot_surface( X_new, Z_new, data_xz_fit, cmap='viridis' )              # 绘制曲面 | Draw the surface
    ax.set_title(  figure_title + "  physical time = " + str(time) )          # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "X [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " xz surface_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " xz surface_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    return

####################################################################################



####################################################################################

## 这是一个对某一二进制数据的画图函数 | This is a plotting function for a given binary data file
## 对 yz 平面数据进行画图 | Plot the data on the yz plane

def get_data_yz( Rmin, Rmax, n, data0, time, figure_title, figure_outdir ):

    figure_contourplot_outdir = os.path.join(figure_outdir, "contour plot")
    figure_densityplot_outdir = os.path.join(figure_outdir, "density plot")
    figure_surfaceplot_outdir = os.path.join(figure_outdir, "surface plot")

    ## 根据读到的格点信息还原格点坐标 | Reconstruct the grid coordinates from the read-in grid information
    x = numpy.linspace(Rmin[0], Rmax[0], n[0])
    y = numpy.linspace(Rmin[1], Rmax[1], n[1])
    z = numpy.linspace(Rmin[2], Rmax[2], n[2])

    ## 用 meshgrid 建立二维格点坐标 | Build the 2D grid coordinates with meshgrid                             
    Y, Z = numpy.meshgrid(y, z)    
    
    ## 补充 numpy.meshgrid 的相关信息 | Supplementary information on numpy.meshgrid
    ## 假设 x 是长度为 nx 的数组，y 为长度为 ny 的数组 | Suppose x is an array of length nx and y is an array of length ny
    ## 而 X, Y = numpy.meshgrid(x, y) 得到的 X 和 Y 都是 (ny, nx) 数组，并且 X 的每一行都是 x 的副本，Y 的每一列都是 y 的副本 | Then X, Y = numpy.meshgrid(x, y) gives X and Y both as (ny, nx) arrays, where each row of X is a copy of x and each column of Y is a copy of y   

    ## 获取 yz 平面上的数据 | Get the data on the yz plane
    '''
    data_yz = data0[:,:,0]
    '''
    index = 0
    ## 找出最接近 x=0 的数据 | Find the data closest to x=0
    for i in range( n[0] ):
        if i > 0:
            if x[i-1] <= 0.0 < x[i]:
                if ( abs(x[i] - 0.0) <= abs(x[i-1] - 0.0) ):
                    index = i
                else:
                    index = i - 1
    data_yz = data0[:,:,index]
                
        
    ## 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    ## 经过测试后发现无需转置 | Tests showed that transposition is not needed
    ## data_yz_new = numpy.transpose(data_yz)
    
    ## 设置新坐标，方便对数据进行插值 | Set up new coordinates to make interpolating the data easier 
    x_new = numpy.linspace(Rmin[0], Rmax[0], int(2.5*n[0]))
    y_new = numpy.linspace(Rmin[1], Rmax[1], int(2.5*n[1]))
    z_new = numpy.linspace(Rmin[2], Rmax[2], int(2.5*n[2]))
    Y_new, Z_new = numpy.meshgrid(y_new, z_new)
    
    ## 对数据进行插值 | Interpolate the data
    data_yz_fit = scipy.interpolate.griddata( (Y.flatten(), Z.flatten()), data_yz.flatten(), (Y_new, Z_new), method="cubic" )

    # 下面画出二维等高线图 | Below, draw the 2D contour plot
    fig, ax = plt.subplots()
    # contourf = ax.contourf( Y, Z, data_yz, 8, cmap='coolwarm', norm=LogNorm(vmin=1, vmax=10), levels=numpy.logspace(-2, 2, 8))  # 使用'coolwarm'色板，并设置标准色彩映射 | Use the 'coolwarm' colormap and set the standard color mapping
    # contourf = ax.contourf( Y, Z, data_yz_0, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  Y, Z, data_yz_0, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    # 由于原始数据是按照列排列的，因此我们转置之后再画图 | Since the raw data are arranged in columns, we transpose them before plotting
    # contourf = ax.contourf( Y, Z, data_yz, cmap=plt.get_cmap('RdYlGn_r') )
    # contour  = ax.contour(  Y, Z, data_yz, 8, colors='k', linewidths=0.5 )       # 添加等高线 | Add contour lines
    # 对插值后的数据进行画图 | Plot the interpolated data
    contourf = ax.contourf( Y_new, Z_new, data_yz_fit, cmap=plt.get_cmap('RdYlGn_r') )
    contour  = ax.contour(  Y_new, Z_new, data_yz_fit, 8, colors='k', linewidths=0.5 )     # 添加等高线 | Add contour lines
    cbar     = plt.colorbar(contourf)                                                      # 添加色条 | Add a colorbar
    ax.set_title(  figure_title + "  physical time = " + str(time) )                       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "Y [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " yz contour_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_contourplot_outdir, figure_title + " time = " + str(time) + " yz contour_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()
    
    # 下面画出二维热图 | Below, draw the 2D heat map    
    # fig1 = plt.figure()
    fig1, ax  = plt.subplots()
    # 经过测试后发现，似乎不用转置，直接用 imshowfig 画图就可以得到相应的 yz 图 | Tests showed that transposition seems unnecessary: plotting directly with imshowfig gives the correct yz figure
    # imshowfig = plt.imshow( data_yz, interpolation='bicubic', extent=[Y.min(), Y.max(), Z.min(), Z.max()] )
    # imshowfig = plt.imshow( numpy.transpose(data_yz), interpolation='bicubic', extent=[Y.min(), Y.max(), Z.min(), Z.max()] )
    # 不知道为什么，z轴坐标是反的，需要人为处理以下才能化成正确的图 | For some reason the z-axis coordinates are reversed, and manual handling is needed to get a correct figure
    imshowfig = plt.imshow( data_yz, interpolation='bicubic', extent=[Y.min(), Y.max(), Z.max(), Z.min()] )                                                    
    cbar      = plt.colorbar(imshowfig)                                     # 添加色条 | Add a colorbar
    ax.invert_yaxis()                                                       # 将 y 轴的正负翻转（对应 z 的数据） | Flip the sign of the y axis (for the z data)
    ax.set_title(  figure_title + "  physical time = " + str(time)  )       # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "Y [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " yz density_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_densityplot_outdir, figure_title + " time = " + str(time) + " yz density_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    # 下面画出三维图 | Below, draw the 3D figure
    fig2 = plt.figure()                                                       # 创建一个新的图像 | Create a new figure
    ax = fig2.add_subplot( 111, projection='3d' )                             # 创建一个 3D 绘图区域 | Create a 3D plotting area
    # 对插值后的数据进行画图 | Plot the interpolated data
    # ax.plot_surface( X, Z, data_xz, cmap='viridis' )                        # 绘制曲面 | Draw the surface
    ax.plot_surface( Y_new, Z_new, data_yz_fit, cmap='viridis' )              # 绘制曲面 | Draw the surface
    ax.set_title(  figure_title + "  physical time = " + str(time) )          # 设置标题和轴标签 | Set the title and axis labels
    ax.set_xlabel( "Y [M]" )
    ax.set_ylabel( "Z [M]" )
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " yz surface_plot.pdf") )   # 保存图像 | Save the figure
    plt.savefig( os.path.join(figure_surfaceplot_outdir, figure_title + " time = " + str(time) + " yz surface_plot.png") )   # 由于后面需要使用cv2制作动画，再保存一份png图像 | Since cv2 is needed later to make animations, also save a png copy
    plt.close()

    return

####################################################################################



####################################################################################

## 该函数生成若干文件夹存放二进制数据对应的图 | This function creates several folders to store the figures of the binary data

def generate_binary_data_plot_directionary( figure_outdir ):

    ## 设定文件夹目录 | Set the folder directories
    figure_outdir_xy = os.path.join( figure_outdir, "XY plot" )
    figure_outdir_xz = os.path.join( figure_outdir, "XZ plot" )
    figure_outdir_yz = os.path.join( figure_outdir, "YZ plot" )
    
    surface_plot_outdir_xy = os.path.join( figure_outdir_xy, "surface plot" )
    surface_plot_outdir_xz = os.path.join( figure_outdir_xz, "surface plot" )
    surface_plot_outdir_yz = os.path.join( figure_outdir_yz, "surface plot" )
    
    density_plot_outdir_xy = os.path.join( figure_outdir_xy, "density plot" )
    density_plot_outdir_xz = os.path.join( figure_outdir_xz, "density plot" )
    density_plot_outdir_yz = os.path.join( figure_outdir_yz, "density plot" )

    contour_plot_outdir_xy = os.path.join( figure_outdir_xy, "contour plot" )
    contour_plot_outdir_xz = os.path.join( figure_outdir_xz, "contour plot" )
    contour_plot_outdir_yz = os.path.join( figure_outdir_yz, "contour plot" )
    
    ## 根据输入文件设定生成对应的文件夹 | Generate the corresponding folders according to the input file settings

    if (input_data.plot_binary_data_set == "xy-xz-yz-plot"):
        
        os.mkdir( figure_outdir_xy )
        os.mkdir( figure_outdir_xz )
        os.mkdir( figure_outdir_yz )

        os.mkdir( surface_plot_outdir_xy )
        os.mkdir( surface_plot_outdir_xz )
        os.mkdir( surface_plot_outdir_yz )

        os.mkdir( density_plot_outdir_xy )
        os.mkdir( density_plot_outdir_xz )
        os.mkdir( density_plot_outdir_yz )

        os.mkdir( contour_plot_outdir_xy )
        os.mkdir( contour_plot_outdir_xz )
        os.mkdir( contour_plot_outdir_yz )

    elif (input_data.plot_binary_data_set == "xy-plot"):

        os.mkdir( figure_outdir_xy )

        os.mkdir( surface_plot_outdir_xy )
        os.mkdir( density_plot_outdir_xy )
        os.mkdir( contour_plot_outdir_xy )

    elif (input_data.plot_binary_data_set == "xz-plot"):

        os.mkdir( figure_outdir_xz )

        os.mkdir( surface_plot_outdir_xz )
        os.mkdir( density_plot_outdir_xz )
        os.mkdir( contour_plot_outdir_yz )

    elif (input_data.plot_binary_data_set == "yz-plot"):

        os.mkdir( figure_outdir_yz )

        os.mkdir( surface_plot_outdir_yz )
        os.mkdir( density_plot_outdir_yz )
        os.mkdir( contour_plot_outdir_yz )

    return

####################################################################################

