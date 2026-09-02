
##################################################################
##
## AMSS-NCKU 数值相对论启动程序
## 小曲
## 2024/03/19
## 2026/02/09 修改
##
##################################################################


##################################################################

## 程序即将开始，选择屏幕终端的输出语言

import logo
import main_program_function

## 输出本程序的 logo
logo.print_logo()

## 让用户选择一语言    
input_language = main_program_function.setting_language() 

##################################################################

## 输出本程序的说明

import print_information

if ( input_language == "Chinese" ):
    print_information.print_program_introduction_Chinese()
elif ( input_language == "English" ):
    print_information.print_program_introduction_English()


## 屏幕输出程序开始运行前的提示
print_information.print_begin_program( input_language )
## 如果同意继续运算，则输入回车
inputvalue = input()   ## 设定一个输入（回车），以便程序下一步运行         
print()


##################################################################

## 导入参数文件

import AMSS_NCKU_Input as input_data


##################################################################

## 生成文件目录，用来储存程序运行的数据

import os
import shutil


## 根据输入文件来设置文件目录
File_directory = os.path.join(input_data.File_directory)   

## 检查设定的输出目录是否存在，并根据用户的输入决定继续计算还是终止计算
main_program_function.continue_or_stop( File_directory, input_language )
        
## 如果文件目录已经存在，则去除它
shutil.rmtree(File_directory, ignore_errors=True)

## 创建文件夹
os.mkdir(File_directory)

## 复制 python 输入文件到程序运行目录
shutil.copy("AMSS_NCKU_Input.py", File_directory)

## 生成文件目录，用来存放各类输出文件

output_directory = os.path.join(File_directory, "AMSS_NCKU_output")
os.mkdir(output_directory)

binary_results_directory = os.path.join(output_directory, input_data.Output_directory)
os.mkdir(binary_results_directory)

figure_directory = os.path.join(File_directory, "figure")
os.mkdir(figure_directory)

video_directory = os.path.join(File_directory, "video")
os.mkdir(video_directory)

## 屏幕输出文件目录已生成
print_information.print_make_directory( input_language )


##################################################################

## 输出相关的参数信息

import setup

setup_data = setup.setup( input_language ) 

## 输出相关的参数信息
if ( input_language == "Chinese" ):
    setup.print_input_data_Chinese( setup_data, File_directory )
elif ( input_language == "English" ):
    setup.print_input_data_English( setup_data, File_directory )

## 输出分辨率的相关信息
## 如果分辨率满足要求，则进行下一步计算
print_information.print_whether_grid_is_satisfied( input_language )
inputvalue = input()  ## 设定一个输入（回车），以便程序下一步运行
print()


##################################################################

## 初始化并输出 puncture 相关的信息

import puncture_initialize

## 初始化 puncture 数据
puncture_data = puncture_initialize.generate_puncture_input_data( input_language )

## 屏幕输出 puncture 的信息
puncture_initialize.print_puncture_information( input_language, puncture_data )
inputvalue = input()  ## 设定一个输入（回车），以便程序下一步运行
print() 


##################################################################

## 记录程序开始时间

import time

start_time = time.time()  # 记录程序开始时间


##################################################################

## 设定初始时刻的 AMR 网格

import numerical_grid        

## 生成一个初始时刻的AMR网格
initial_AMR_grid_data = numerical_grid.generate_initial_AMR_grid( input_language, puncture_data )

## 画出初始的格点设置
numerical_grid.plot_initial_grid( input_language, initial_AMR_grid_data )


##################################################################

## 自动生成 AMSS-NCKU C++ 程序的输入文件

## 根据设定的参数生成 AMSS-NCKU 程序的输入文件   
setup.generate_AMSSNCKU_input( input_language )

## 根据格点信息生成 cgh 的相关输入
numerical_grid.append_AMSSNCKU_cgh_input( input_language, initial_AMR_grid_data )  

##  根据算法和参数生成 AMSS-NCKU 的宏文件
main_program_function.generate_macrodef_file( input_language )


##################################################################

## 根据用户的要求编译 AMSS-NCKU 程序

## 屏幕输出正在编译 AMSS-NCKU
print_information.print_compile_AMSS_NCKU( input_language )

AMSS_NCKU_source_path  = "AMSS_NCKU_source"
AMSS_NCKU_compile_path = os.path.join(File_directory, "AMSS_NCKU_compile")

###############################

## 如果 AMSS-NCKU 的源文件夹不存在，则创建它

# if not os.path.exists(destination_folder):
#     os.makedirs(destination_folder)

if not os.path.exists(AMSS_NCKU_source_path):
    os.makedirs(AMSS_NCKU_source_path)
    print_information.print_compile_AMSS_NCKU_debug( input_language )  ## 输出编译出错的信息
    
###############################

# 拷贝 AMSS-NCKU 的源文件，准备编译
shutil.copytree(AMSS_NCKU_source_path, AMSS_NCKU_compile_path)    

# 将整个src文件夹复制到dst位置
# shutil.copytree(src, dst)

# 将生成的宏文件拷贝进 AMSS-NCKU 的代码文件夹

macrodef_h_path  = os.path.join(File_directory, "macrodef.h") 
macrodef_fh_path = os.path.join(File_directory, "macrodef.fh") 

shutil.copy2(macrodef_h_path,  AMSS_NCKU_compile_path)
shutil.copy2(macrodef_fh_path, AMSS_NCKU_compile_path)

# 复制文件到目标位置
# shutil.copy2(source_file_path, target_location)
# shutil.copy2保留文件的元数据，如修改时间。如果你只想复制文件内容，不保留元数据，可以使用shutil.copy。

###############################

# 编译相关程序

import makefile_and_run_ZJU as makefile_and_run

## 先切换到目标文件夹
os.chdir(AMSS_NCKU_compile_path)
 
## 编译 AMSS-NCKU 的主程序 ABE/ABEGPU
makefile_and_run.makefile_ABE( input_language )

## 如果初值方法选为 TwoPuncture，编译可执行文件 TwoPunctureABE 
## 如果初值方法选为其它方法，则跳过这一步

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ): 
    makefile_and_run.makefile_TwoPunctureABE( input_language )
else:
    print()
    
## 如果设置 puncture 的方式选为 Automatically-BBH，且允许后牛顿演化，编译可执行文件 PNOrbit 
## 如果设置 puncture 的方式选为 Manually，或 Allow_PN_Evaluation 选为 "no"，则跳过这一步

if ( input_data.puncture_data_set == "Automatically-BBH" and input_data.Allow_PN_Evaluation == "yes" ): 
    makefile_and_run.makefile_PNOrbit( input_language )
else: 
    print()
    
###########################

## 改变当前工作目录到上一级目录
os.chdir('..')
os.chdir('..')

print()

##################################################################

## 将 AMSS-NCKU 程序可执行文件 ABE 复制到程序运行的文件夹

if (input_data.GPU_Calculation == "no"):
    ABE_file = os.path.join(AMSS_NCKU_compile_path, "ABE")
elif (input_data.GPU_Calculation == "yes"):
    ABE_file = os.path.join(AMSS_NCKU_compile_path, "ABEGPU")

if not os.path.exists( ABE_file ):
    print_information.print_no_ABE_error( input_language )  ## 屏幕输出不存在 ABE/ABEGPU 的报错信息

## 复制可执行文件 ABE 到程序运行目录
shutil.copy2(ABE_file, output_directory)

###########################

## 如果初值方法选为 TwoPuncture，将 AMSS-NCKU 程序可执行文件 TwoPunctureABE 复制到程序运行的文件夹

TwoPuncture_file = os.path.join(AMSS_NCKU_compile_path, "TwoPunctureABE")

if ( input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):

    ## 如果找不到可执行文件 TwoPunctureABE，屏幕输出报错信息
    if not os.path.exists( TwoPuncture_file ):
        print_information.print_no_TwoPunture_error( input_language ) 

    ## 复制可执行文件 TwoPunctureABE 到程序运行目录
    shutil.copy2(TwoPuncture_file, output_directory)

###########################

## 如果设置 puncture 的方式选为 Automatically-BBH，且允许后牛顿演化
## 将 AMSS-NCKU 程序可执行文件 PNOrbit 拷贝到程序运行的文件夹

PNOrbit_file = os.path.join(AMSS_NCKU_compile_path, "PNOrbit")

if ( input_data.puncture_data_set == "Automatically-BBH" and input_data.Allow_PN_Evaluation == "yes" ):
    
    ## 如果找不到可执行文件 PNOrbit，屏幕输出报错信息
    if not os.path.exists( PNOrbit_file ):
        print_information.print_no_PNOrbit_error( input_language ) 

    ## 复制可执行文件 PNOrbit 到程序运行目录
    shutil.copy2(PNOrbit_file, output_directory)

##################################################################

## 如果设置 puncture 的方式选为 Automatically-BBH，且允许后牛顿演化，生成并运行 PNOrbit 程序
## 用后牛顿演化计算更精确的双星轨道参数（从 Distance_initial 演化到 Distance_final）
## 如果 Allow_PN_Evaluation 选为 "no"，则跳过 PNOrbit

if ( input_data.puncture_data_set == "Automatically-BBH" and input_data.Allow_PN_Evaluation == "yes" ):

    ## 屏幕输出生成 PNOrbit 输入文件
    print_information.print_generate_PNOrbit_input( input_language )
    
    import generate_PNOrbit_input
    
    ## 生成 PNOrbit 的输入文件
    generate_PNOrbit_input.generate_AMSSNCKU_PNOrbit_input( input_language, puncture_data )
    print_information.print_finish_PNOrbit_input( input_language )   ## 屏幕输出 PNOrbit 的输入文件已完成
    
    ## 生成的 PNOrbit 输入文件名
    AMSS_NCKU_PNOrbit_inputfile      = 'AMSS-NCKU-PNOrbit.input'
    AMSS_NCKU_PNOrbit_inputfile_path = os.path.join( File_directory, AMSS_NCKU_PNOrbit_inputfile )
 
    ## 复制并重命名文件（PNOrbit 程序读取当前目录下的 PN_Orbit.par）
    shutil.copy2( AMSS_NCKU_PNOrbit_inputfile_path, os.path.join(output_directory, 'PN_Orbit.par') )
    
    ###########################

    ## 运行 PNOrbit 来生成更精确的轨道参数
    
    print()
    ## print( " 准备启动 AMSS-NCKU PNOrbit 程序 ，按回车继续！！！ " )
    ## inputvalue = input()                    
    print()
    
    ## 先切换到目标文件夹
    os.chdir(output_directory)

    ## 运行 PNOrbit 程序
    makefile_and_run.run_PNOrbit( input_language )
    
    ###########################
    
    ## 改变当前工作目录到上两级目录
    os.chdir('..')
    os.chdir('..')
    
    ###########################

    ## 根据 PNOrbit 的运行结果更新 puncture data（解析 PNorbit.dat 并更新位置、动量、间距）
    ## import puncture_initialize
    
    puncture_data_new = puncture_initialize.update_puncture_data_by_PNOrbit( input_language, puncture_data )

else:

    puncture_data_new = puncture_data
    
##################################################################

## 如果初值方法选为 TwoPuncture，生成 TwoPuncture 的输入文件

if ( input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
    
    ## 屏幕输出生成 TwoPuncture 输入文件
    print_information.print_generate_TwoPunture_input( input_language )
    
    import generate_TwoPuncture_input
    
    ## 根据双星模型设定 TwoPuncture 需要的数据
    ## 新版已放入 puncture_initilize 中
    ## puncture_data_new = generate_TwoPuncture_input.generate_puncture_input_data( input_language )
    
    ## 生成 TwoPuncture 的输入文件
    generate_TwoPuncture_input.generate_AMSSNCKU_TwoPuncture_input( input_language, puncture_data_new )
    print_information.print_finish_TwoPunture_input( input_language )   ## 屏幕输出 TwoPuncture 的输入文件已完成
    
    ## 生成的 AMSS-NCKU 输入文件名
    AMSS_NCKU_TwoPuncture_inputfile      = 'AMSS-NCKU-TwoPuncture.input'
    AMSS_NCKU_TwoPuncture_inputfile_path = os.path.join( File_directory, AMSS_NCKU_TwoPuncture_inputfile )
 
    ## 复制并重命名文件
    shutil.copy2( AMSS_NCKU_TwoPuncture_inputfile_path, os.path.join(output_directory, 'TwoPunctureinput.par') )
    
    ###########################

    ## 运行 TwoPuncture 来生成初值文件
    
    print()
    ## print( " 准备启动 AMSS-NCKU TwoPuncture 程序 ，按回车继续！！！ " )
    ## inputvalue = input()                    
    print()
    
    ## 先切换到目标文件夹
    os.chdir(output_directory)

    ## 运行 TwoPuncture 程序
    makefile_and_run.run_TwoPunctureABE( input_language )
    
    ###########################
    
    ## 改变当前工作目录到上两级目录
    os.chdir('..')
    os.chdir('..')
    
##################################################################
    
## 根据 TwoPuncture 的运行结果来更新 puncture data

import renew_puncture_parameter
    
renew_puncture_parameter.append_AMSSNCKU_BSSN_input( puncture_data_new, File_directory, output_directory )


## 生成的 AMSS-NCKU 输入文件名
AMSS_NCKU_inputfile      = 'AMSS-NCKU.input'
AMSS_NCKU_inputfile_path = os.path.join(File_directory, AMSS_NCKU_inputfile)
 
## 复制并重命名文件
shutil.copy2( AMSS_NCKU_inputfile_path, os.path.join(output_directory, 'input.par') )

## 屏幕输出成功拷贝参数文件
print_information.copy_input_parfile( input_language )
    

##################################################################

## 准备启动 AMSS-NCKU 的可执行程序 ABE/ABEGPU
print_information.print_begin_ABE( input_language )

## 先切换到目标文件夹
os.chdir( output_directory )

## 执行 ABE/ABEGPU
makefile_and_run.run_ABE( input_language )

## 改变当前工作目录到上两级目录
os.chdir('..')
os.chdir('..')


##################################################################

## 将一些基本输入拷贝出来，方便进行 debug

## 用于输出计算所用参数的文件地址
AMSS_NCKU_error_file_path = os.path.join(binary_results_directory, "setting.par")
## 复制并重命名文件
shutil.copy( AMSS_NCKU_error_file_path, os.path.join(output_directory, "AMSSNCKU_setting_parameter") )

## 用于报错的文件地址
AMSS_NCKU_error_file_path = os.path.join(binary_results_directory, "Error.log")
## 复制并重命名文件
shutil.copy( AMSS_NCKU_error_file_path, os.path.join(output_directory, "Error.log") )

## 程序基本输出
AMSS_NCKU_BH_data         = os.path.join(binary_results_directory, "bssn_BH.dat"        )
AMSS_NCKU_ADM_data        = os.path.join(binary_results_directory, "bssn_ADMQs.dat"     )
AMSS_NCKU_psi4_data       = os.path.join(binary_results_directory, "bssn_psi4.dat"      )
AMSS_NCKU_constraint_data = os.path.join(binary_results_directory, "bssn_constraint.dat")
## 复制并重命名文件
shutil.copy( AMSS_NCKU_BH_data,         os.path.join(output_directory, "bssn_BH.dat"        ) )
shutil.copy( AMSS_NCKU_ADM_data,        os.path.join(output_directory, "bssn_ADMQs.dat"     ) )
shutil.copy( AMSS_NCKU_psi4_data,       os.path.join(output_directory, "bssn_psi4.dat"      ) )
shutil.copy( AMSS_NCKU_constraint_data, os.path.join(output_directory, "bssn_constraint.dat") )

## 程序其它输出
if (input_data.Equation_Class == "BSSN-EM"):
    AMSS_NCKU_phi1_data = os.path.join(binary_results_directory, "bssn_phi1.dat" )
    AMSS_NCKU_phi2_data = os.path.join(binary_results_directory, "bssn_phi2.dat" )
    shutil.copy( AMSS_NCKU_phi1_data, os.path.join(output_directory, "bssn_phi1.dat" ) )
    shutil.copy( AMSS_NCKU_phi2_data, os.path.join(output_directory, "bssn_phi2.dat" ) )
elif (input_data.Equation_Class == "BSSN-EScalar"):
    AMSS_NCKU_maxs_data = os.path.join(binary_results_directory, "bssn_maxs.dat" )
    shutil.copy( AMSS_NCKU_maxs_data, os.path.join(output_directory, "bssn_maxs.dat" ) )

##################################################################

##  对 AMSS-NCKU 程序运行结果进行画图

import plot_xiaoqu
import generate_video_xiaoqu
import plot_GW_strain_amplitude_xiaoqu

## 屏幕输出开始画图
print_information.print_begin_plot( input_language )

## 画出黑洞轨迹图
plot_xiaoqu.generate_puncture_orbit_plot(   input_language, binary_results_directory, figure_directory )
plot_xiaoqu.generate_puncture_orbit_plot3D( input_language, binary_results_directory, figure_directory )

## 画出黑洞间距随时间的变化图
plot_xiaoqu.generate_puncture_distence_plot( input_language, binary_results_directory, figure_directory )

## 画出引力波波形图
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_gravitational_wave_psi4_plot( input_language,binary_results_directory, figure_directory, i )
    plot_GW_strain_amplitude_xiaoqu.generate_gravitational_wave_amplitude_plot( input_language, binary_results_directory, figure_directory, i )

## 画出时空 ADM 质量的变化图
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_ADMmass_plot( input_language, binary_results_directory, figure_directory, i )

## 画出哈密顿约束违反性况的变化图
for i in range(input_data.grid_level):
    plot_xiaoqu.generate_constraint_check_plot( input_language, binary_results_directory, figure_directory, i )

## 对储存的二进制数据画出图像
plot_xiaoqu.generate_binary_data_plot( input_language, binary_results_directory, figure_directory )

## 对储存的二进制数据生成视频
generate_video_xiaoqu.generate_binary_data_video( input_language, figure_directory, video_directory )

##################################################################

## 输出程序运行时间

end_time = time.time()                # 记录结束时间
elapsed_time = end_time - start_time  # 计算运行时间

## 屏幕输出程序运行时间
print_information.print_time_cost( elapsed_time, input_language )


##################################################################

## 屏幕输出程序结束
print_information.print_program_end( input_language )

##################################################################


