
#################################################
##
## 这个文件用于生成 AMSS-NCKU 数值相对论计算结果的视频
## 小曲
## 2026/02/27
##
#################################################

import os
import shutil
import gc
import cv2
import numpy
import re

import AMSS_NCKU_Input as input_data



## 图片帧率（每秒显示的图片数量）
## 你可以根据需要调整这个值

fps_set = 4


####################################################################################

## 该函数根据二进制数据画出所有二维图

def generate_binary_data_video( input_language, figure_outdir, video_outdir ):


    if ( input_language == "Chinese" ):
        print(                                             )
        print( " 开始对 AMSS-NCKU 程序输出的二进制数据制作动画 " )
        print(                                             )
    elif ( input_language == "English" ):
        print(                                                               )
        print( " Beginning the AMSS-NCKU Binary Data Plotting From Figures " )
        print(                                                               )

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

    ## 设定产生视频的文件夹目录
    
    ## 利用 numpy 设置空数组，用于存放字符串
    ## 设置 dtype='U200'，确保存放的是字符串，长度为 200（防止长度不够）
    video_outdir_static  = numpy.empty( input_data.static_grid_level, dtype='U200' )
    video_outdir_moving  = numpy.empty( input_data.moving_grid_level, dtype='U200' )
    video_outdir_moving2 = numpy.empty( (input_data.moving_grid_level, input_data.puncture_number), dtype='U200' )
    
    ## 设置每层网格画图文件夹的目录
    for i in range(input_data.static_grid_level):
        video_outdir_static[i] = os.path.join( video_outdir, "level" + str(i) )
    for i in range(input_data.moving_grid_level):
        j = i + input_data.static_grid_level
        for k in range(input_data.puncture_number):
            video_outdir_moving[i] = os.path.join( video_outdir, "level" + str(j) )
            video_outdir_moving2[i,k] = os.path.join( video_outdir_moving[i], "puncture" + str(k) )


    ###########################################

    ## 如果输入文件中 plot_binary_data_level 设定为 "All-Level"，对所有网格层的二进制和数据制作动画

    if (input_data.plot_binary_data_level == "All-Level"): 
        
        ## 对固定网格层的二进制数制作动画
        for i in range(input_data.static_grid_level):
            
            ## 生成该层二进制数据制作动画的文件夹目录
            os.mkdir( video_outdir_static[i] )
            
            ## 生成该层二进制数据制作动画的各种子文件夹目录
            generate_binary_data_video_directionary( video_outdir_static[i] )

            ## 对该固定网格层的二进制数据制作动画
            generate_vedio_from_each_level( input_language, i, 0, figure_outdir_static[i], video_outdir_static[i] )

        ## 对移动网格层的二进制数制作动画
        for i in range(input_data.moving_grid_level):

            ## 生成该层二进制数据制作动画的文件夹目录
            os.mkdir( video_outdir_moving[i] )
            
            j = i + input_data.static_grid_level
            
            for k in range(input_data.puncture_number):

                ## 生成该层二进制数据制作动画的各种子文件夹目录
                os.mkdir( video_outdir_moving2[i,k] )
                generate_binary_data_video_directionary( video_outdir_moving2[i,k] )

                ## 对该移动网格层的二进制数据制作动画
                generate_vedio_from_each_level( input_language, i, k, figure_outdir_moving2[i,k], video_outdir_moving2[i,k] )

    ###########################################

    ## 如果输入文件中 plot_binary_data_level 设定为 "Single-Level"，对设定的网格层的二进制和数据制作动画

    elif (input_data.plot_binary_data_level == "Single-Level"):

        ## 设定的网格层在固定网格层
        if (input_data.plot_binary_data_levelnumber < input_data.static_grid_level):

            ## 生成该层二进制数据制作动画的文件夹目录
            os.mkdir( video_outdir_static[input_data.plot_binary_data_levelnumber] )
            
            ## 生成该层二进制数据制作动画的各种子文件夹目录
            generate_binary_data_video_directionary( video_outdir_static[i] )

            ## 对该固定网格层的二进制数据制作动画
            generate_vedio_from_each_level( input_language, input_data.plot_binary_data_levelnumber, figure_outdir_static[input_data.plot_binary_data_levelnumber] )

        ## 设定的网格层在固定网格层
        else:
            
            j = input_data.plot_binary_data_levelnumber - input_data.static_grid_level

            ## 生成该层二进制数据制作动画的文件夹目录
            os.mkdir( video_outdir_moving[j] )
            
            for k in range(input_data.puncture_number):

                ## 生成该层二进制数据制作动画的各种子文件夹目录
                os.mkdir( video_outdir_moving2[j,k] )
                generate_binary_data_video_directionary( video_outdir_moving2[j,k] )

                ## 对该移动网格层的二进制数据制作动画
                generate_vedio_from_each_level( input_language, j, k, figure_outdir_moving2[j,k], video_outdir_moving2[j,k] )


    ###########################################

    if ( input_language == "Chinese" ):
        print(                          )
        print( " 二进制数据制作动画已完成 " )
        print(                          )
    elif ( input_language == "English" ):
        print(                                              )
        print( " Video for Binary Data Has been Generated " )
        print(                                              )

    return

####################################################################################




####################################################################################

## 该函数生成若干文件夹存放二进制数据对应的动画

def generate_binary_data_video_directionary( video_outdir ):

    ## 设定文件夹目录
    video_outdir_xy = os.path.join( video_outdir, "XY video" )
    video_outdir_xz = os.path.join( video_outdir, "XZ video" )
    video_outdir_yz = os.path.join( video_outdir, "YZ video" )
    
    surface_plot_video_outdir_xy = os.path.join( video_outdir_xy, "surface plot video" )
    surface_plot_video_outdir_xz = os.path.join( video_outdir_xz, "surface plot video" )
    surface_plot_video_outdir_yz = os.path.join( video_outdir_yz, "surface plot video" )
    
    density_plot_video_outdir_xy = os.path.join( video_outdir_xy, "density plot video" )
    density_plot_video_outdir_xz = os.path.join( video_outdir_xz, "density plot video" )
    density_plot_video_outdir_yz = os.path.join( video_outdir_yz, "density plot video" )

    contour_plot_video_outdir_xy = os.path.join( video_outdir_xy, "contour plot video" )
    contour_plot_video_outdir_xz = os.path.join( video_outdir_xz, "contour plot video" )
    contour_plot_video_outdir_yz = os.path.join( video_outdir_yz, "contour plot video" )
    
    ## 根据输入文件设定生成对应的文件夹

    if (input_data.plot_binary_data_set == "xy-xz-yz-plot"):
        
        os.mkdir( video_outdir_xy )
        os.mkdir( video_outdir_xz )
        os.mkdir( video_outdir_yz )

        os.mkdir( surface_plot_video_outdir_xy )
        os.mkdir( surface_plot_video_outdir_xz )
        os.mkdir( surface_plot_video_outdir_yz )

        os.mkdir( density_plot_video_outdir_xy )
        os.mkdir( density_plot_video_outdir_xz )
        os.mkdir( density_plot_video_outdir_yz )

        os.mkdir( contour_plot_video_outdir_xy )
        os.mkdir( contour_plot_video_outdir_xz )
        os.mkdir( contour_plot_video_outdir_yz )

    elif (input_data.plot_binary_data_set == "xy-plot"):

        os.mkdir( video_outdir_xy )

        os.mkdir( surface_plot_video_outdir_xy )
        os.mkdir( density_plot_video_outdir_xy )
        os.mkdir( contour_plot_video_outdir_xy )

    elif (input_data.plot_binary_data_set == "xz-plot"):

        os.mkdir( video_outdir_xz )

        os.mkdir( surface_plot_video_outdir_xz )
        os.mkdir( density_plot_video_outdir_xz )
        os.mkdir( contour_plot_video_outdir_yz )

    elif (input_data.plot_binary_data_set == "yz-plot"):

        os.mkdir( video_outdir_yz )

        os.mkdir( surface_plot_video_outdir_yz )
        os.mkdir( density_plot_video_outdir_yz )
        os.mkdir( contour_plot_video_outdir_yz )

    return

####################################################################################




####################################################################################

## 这个函数对某一个网格层中的数据制作动画
## 函数输入为
## 网格层编号      level_number
## puncture 编号  puncture_number
## 该网格层对应的图片文件夹 figure_outdir
## 该网格层对应的视频文件夹 video_outdir

def generate_vedio_from_each_level( input_language, level_number, puncture_number, figure_outdir, video_outdir ):

    
    ## 屏幕输出开始生成该网格层的视频
    if input_language == "Chinese":
        if level_number < input_data.static_grid_level:
            print( " 正在生成固定网格层二进制数据对应的视频 level number =", level_number )
            print(                                                                 )
        else:
            print( " 正在生成移动网格层二进制数据对应的视频 level number =", level_number, " puncture number =", puncture_number )
            print(                                                                                                        )
    elif input_language == "English":
        if level_number < input_data.static_grid_level:
            print( " Generating Videos for Binary Data in Static Level, level number =", level_number )
            print(                                                                                    )
        else:
            print( " Generating Videos for Binary Data in Moving Level, level number =", level_number, " puncture number =", puncture_number )
            print(                                                                                                                           )

 
    ## 设定图片文件夹目录
    figure_outdir_xy = os.path.join( figure_outdir, "XY plot" )
    figure_outdir_xz = os.path.join( figure_outdir, "XZ plot" )
    figure_outdir_yz = os.path.join( figure_outdir, "YZ plot" )

    surface_plot_figure_outdir_xy = os.path.join( figure_outdir_xy, "surface plot" )
    surface_plot_figure_outdir_xz = os.path.join( figure_outdir_xz, "surface plot" )
    surface_plot_figure_outdir_yz = os.path.join( figure_outdir_yz, "surface plot" )
    
    density_plot_figure_outdir_xy = os.path.join( figure_outdir_xy, "density plot" )
    density_plot_figure_outdir_xz = os.path.join( figure_outdir_xz, "density plot" )
    density_plot_figure_outdir_yz = os.path.join( figure_outdir_yz, "density plot" )

    contour_plot_figure_outdir_xy = os.path.join( figure_outdir_xy, "contour plot" )
    contour_plot_figure_outdir_xz = os.path.join( figure_outdir_xz, "contour plot" )
    contour_plot_figure_outdir_yz = os.path.join( figure_outdir_yz, "contour plot" )

    ## 设定视频文件夹目录
    video_outdir_xy = os.path.join( video_outdir, "XY video" )
    video_outdir_xz = os.path.join( video_outdir, "XZ video" )
    video_outdir_yz = os.path.join( video_outdir, "YZ video" )

    surface_plot_video_outdir_xy = os.path.join( video_outdir_xy, "surface plot video" )
    surface_plot_video_outdir_xz = os.path.join( video_outdir_xz, "surface plot video" )
    surface_plot_video_outdir_yz = os.path.join( video_outdir_yz, "surface plot video" )
    
    density_plot_video_outdir_xy = os.path.join( video_outdir_xy, "density plot video" )
    density_plot_video_outdir_xz = os.path.join( video_outdir_xz, "density plot video" )
    density_plot_video_outdir_yz = os.path.join( video_outdir_yz, "density plot video" )

    contour_plot_video_outdir_xy = os.path.join( video_outdir_xy, "contour plot video" )
    contour_plot_video_outdir_xz = os.path.join( video_outdir_xz, "contour plot video" )
    contour_plot_video_outdir_yz = os.path.join( video_outdir_yz, "contour plot video" )
    
    ## 根据输入文件设定对哪个平面的数据进行画图
    ## 对 xy yz zx 平面数据都画图
    if (input_data.plot_binary_data_set == "xy-xz-yz-plot"):
        ## xy 部分
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_xy, surface_plot_video_outdir_xy )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_xy, density_plot_video_outdir_xy )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_xy, contour_plot_video_outdir_xy )
        ## xz 部分
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_xz, surface_plot_video_outdir_xz )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_xz, density_plot_video_outdir_xz )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_xz, contour_plot_video_outdir_xz )
        ## yz 部分
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_yz, surface_plot_video_outdir_yz )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_yz, density_plot_video_outdir_yz )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_yz, contour_plot_video_outdir_yz )
    ## 只对 xy 平面数据都画图
    elif (input_data.plot_binary_data_set == "xy-plot"):
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_xy, surface_plot_video_outdir_xy )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_xy, density_plot_video_outdir_xy )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_xy, contour_plot_video_outdir_xy )
    ## 只对 xz 平面数据都画图
    elif (input_data.plot_binary_data_set == "xz-plot"):
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_xz, surface_plot_video_outdir_xz )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_xz, density_plot_video_outdir_xz )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_xz, contour_plot_video_outdir_xz )
    ## 只对 yz 平面数据都画图
    elif (input_data.plot_binary_data_set == "yz-plot"):
        generate_video_from_each_folder( input_language, surface_plot_figure_outdir_yz, surface_plot_video_outdir_yz )
        generate_video_from_each_folder( input_language, density_plot_figure_outdir_yz, density_plot_video_outdir_yz )
        generate_video_from_each_folder( input_language, contour_plot_figure_outdir_yz, contour_plot_video_outdir_yz )


    ## 屏幕输出该网格层的视频已经完成
    if input_language == "Chinese":
        if level_number < input_data.static_grid_level:
            print( " 已经生成固定网格层二进制数据对应的视频 level number =", level_number )
            print(                                                                 )
        else:
            print( " 正在生成移动网格层二进制数据对应的视频 level number =", level_number, " puncture number =", puncture_number )
            print(                                                                                                        )
    elif input_language == "English":
        if level_number < input_data.static_grid_level:
            print( " Having Generated Videos for Binary Data in Static Level, level number =", level_number )
            print(                                                                                          )
        else:
            print( " Having Generated Videos for Binary Data in Moving Level, level number =", level_number, " puncture number =", puncture_number )
            print(                                                                                                                                 )
 
    return

####################################################################################



####################################################################################

## 对某各个变量的数据制作动画
## 函数输入
## 图片文件夹 image_folder
## 视频文件夹 video_folder

####################################################################################


def generate_video_from_each_folder( input_language, image_folder, video_folder ):

    ###########################################
    
    ## 设定生成的视频文件名称

    video_name_TrK      = "TrK.mp4"
    video_name_Lapse    = "Lapse_Alpha.mp4"
    video_name_Phi      = "Conformal_Phi.mp4"
    video_name_Hamilton = "Constraint_H.mp4"

    video_output_TrK      = os.path.join( video_folder, "TrK.mp4"           )
    video_output_Lapse    = os.path.join( video_folder, "Lapse_Alpha.mp4"   )
    video_output_Phi      = os.path.join( video_folder, "Conformal_Phi.mp4" )
    video_output_Hamilton = os.path.join( video_folder, "Constraint_H.mp4"  )

    ###########################################

    # 获取图片列表
    image_files =  [ f for f in os.listdir(image_folder) if f.endswith(('.png', '.jpg', '.jpeg')) ]

    ## 筛选出关于各个变量的数据图片
    image_files_TrK      =  [ f for f in image_files if f.startswith(('trK')) ]
    image_files_Lapse    =  [ f for f in image_files if f.startswith(('Lap')) ]
    image_files_Phi      =  [ f for f in image_files if f.startswith(('phi')) ]
    image_files_Hamilton =  [ f for f in image_files if f.startswith(('Cons_Ham')) ]

    ###########################################

    ## 对数据图片进行排序，这样生成视频中才是沿着时间前进方向
    ## 可利用文件名中包含 time = xx.xx 来进行排序

    ## 对数据图片进行排序
    ## image_files_TrK.sort( key=lambda x: int(x) ) 表示对 x 进行排序，int(x) 将 x转换为整数
    ## 但由于图片文件名中又有文本又有数字，所以要利用re.search来提取数字部分

    ## 以下是一个例子
    '''
    text = "aabbss23.226a.pdf"
    # 使用正则表达式查找所有数字，包括小数点
    match = re.search(r'\d+\.\d+', text)
    '''

    image_files_TrK.sort(      key=lambda x: float( re.search(r'\d+\.\d+', x[6:]).group() ) )  ## 为了防止前面的 trk0 进入，从第6位开始索引
    image_files_Lapse.sort(    key=lambda x: float( re.search(r'\d+\.\d+', x[6:]).group() ) )  ## 为了防止前面的 Lap0 进入，从第6位开始索引
    image_files_Phi.sort(      key=lambda x: float( re.search(r'\d+\.\d+', x[6:]).group() ) )  ## 为了防止前面的 phi0 进入，从第6位开始索引
    image_files_Hamilton.sort( key=lambda x: float( re.search(r'\d+\.\d+', x[6:]).group() ) )

    ## 只使用 re.search(r'\d+\.\d+', x[6:]) 是不行的
    ## TypeError: float() argument must be a string or a real number, not 're.Match'
    ## 我们使用 re.search() 来查找第一个匹配项。如果找到了匹配项，我们使用 match.group() 来获取匹配的字符串，然后将其转换为浮点数

    ###########################################

    ## 对 TrK 的数据进行制作动画

    # 获取第一张图片的尺寸，所有图片都应该具有相同的尺寸
    frame                 = cv2.imread( os.path.join(image_folder, image_files_TrK[0]) )
    height, width, layers = frame.shape
    size                  = (width, height)

    # 创建VideoWriter对象
    out = cv2.VideoWriter( video_output_TrK, cv2.VideoWriter_fourcc(*'mp4v'), fps_set, size )

    ## 将图片写入每一帧
    for image_file in image_files_TrK:
        img_path = os.path.join( image_folder, image_file )
        frame    = cv2.imread( img_path )
        out.write( frame )     # 写入每一帧

    out.release()   # 释放资源

    if input_language == "Chinese":
        print( " 视频已生成: ", video_output_TrK )
    elif input_language == "English":
        print( " Video has been generated: ", video_output_TrK )

    ###########################################

    ## 对 Lapse 函数 alpha 的数据进行制作动画

    ## 获取第一张图片的尺寸，所有图片都应该具有相同的尺寸
    frame                 = cv2.imread( os.path.join(image_folder, image_files_Lapse[0]) )
    height, width, layers = frame.shape
    size                  = (width, height)

    # 创建VideoWriter对象
    out = cv2.VideoWriter( video_output_Lapse, cv2.VideoWriter_fourcc(*'mp4v'), fps_set, size )

     ## 将图片写入每一帧
    for image_file in image_files_Lapse:
        img_path = os.path.join( image_folder, image_file )
        frame    = cv2.imread( img_path )
        out.write( frame )    # 写入帧

    out.release()   # 释放资源

    if input_language == "Chinese":
        print( " 视频已生成: ", video_name_Lapse )
    elif input_language == "English":
        print( " Video has been generated: ", video_name_Lapse )

    ###########################################

    ## 对共形因子 phi 的数据进行制作动画

    ## 获取第一张图片的尺寸，所有图片都应该具有相同的尺寸
    frame                 = cv2.imread( os.path.join(image_folder, image_files_Phi[0]) )
    height, width, layers = frame.shape
    size                  = (width, height)

    # 创建VideoWriter对象
    out = cv2.VideoWriter( video_output_Phi, cv2.VideoWriter_fourcc(*'mp4v'), fps_set, size )

     ## 将图片写入每一帧
    for image_file in image_files_Phi:
        img_path = os.path.join( image_folder, image_file )
        frame    = cv2.imread( img_path )
        out.write( frame )    # 写入帧

    out.release()   # 释放资源

    if input_language == "Chinese":
        print( " 视频已生成: ", video_name_Phi )
    elif input_language == "English":
        print( " Video has been generated: ", video_name_Phi )

    ###########################################

    ## 对哈密顿约束的数据进行制作动画

    ## 获取第一张图片的尺寸，所有图片都应该具有相同的尺寸
    frame                 = cv2.imread( os.path.join(image_folder, image_files_Hamilton[0]) )
    height, width, layers = frame.shape
    size                  = (width, height)

    # 创建VideoWriter对象
    out = cv2.VideoWriter( video_output_Hamilton, cv2.VideoWriter_fourcc(*'mp4v'), fps_set, size )

     ## 将图片写入每一帧
    for image_file in image_files_Hamilton:
        img_path = os.path.join( image_folder, image_file )
        frame    = cv2.imread( img_path )
        out.write( frame )    # 写入帧

    out.release()   # 释放资源

    if input_language == "Chinese":
        print( " 视频已生成: ", video_name_Hamilton )
    elif input_language == "English":
        print( " Video has been generated: ", video_name_Hamilton )

    ###########################################

    return

####################################################################################


####################################################################################

# 单独使用的例子
'''
outdir = "./BBH_test"

figure_directory = os.path.join(outdir, "figure")
video_directory  = os.path.join(outdir, "video")

os.mkdir( video_directory )

generate_binary_data_video( "Chinese", figure_directory, video_directory )
'''
####################################################################################
