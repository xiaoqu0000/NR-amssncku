
##################################################################
##
## AMSS-NCKU 数值相对论启动程序
## 小曲
## 2024/03/19
## 2025/12/09 修改
##
##################################################################


##################################################################

## 输出本程序的说明

import print_information

print_information.print_program_introduction()

##################################################################

## 程序开始运行前的提示

print(                                                                                )
print( " 计算即将开始，请确认在 AMSS_NCKU_Input.py 中设置了正确的参数，按回车继续！！！  "     )
print( " 如果输入参数没有设置好，Ctrl+C 退出，调整 AMSS_NCKU_Input.py 中的输入参数！！！ "    )
print(                                                                                )
print( " Simulation will be started, please confirm you have set the correct parameters in the script file " )
print( " AMSS_NCKU_Input.py "                                                                                )
print( " If parameters have been set correctly, press Enter to continue !!!  "                               )
print( " If you have not set parameters，press Ctrl+C to abort the simulation and adjust the parameters "    )
print( " in script file AMSS_NCKU_Input.py !!! "                                                             )
     
## 设定一个输入（回车），以便程序下一步运行
inputvalue = input()           
print()


##################################################################

## 导入参数文件

import AMSS_NCKU_Input as input_data


##################################################################

## 生成文件目录，用来储存程序运行的数据

import os
import shutil
import sys
import time

## 根据输入文件来设置文件目录
File_directory = os.path.join(input_data.File_directory)   

## 如果设定的输出目录存在，则根据用户的选择看是否继续计算
if os.path.exists(File_directory):
    print(                                                                                           )
    print( " 设定的输出目录存在，是否同意覆盖当前目录？ "                                                  )
    print( " 如果同意覆当前目录，请输入'continue'继续计算 "                                               )
    print( " 如果不同意覆当前目录，请输入'stop'退出计算，并在输入文件 AMSS_NCKU_Input.py 中重新设置输出目录 "  )
    print( " Output dictionary has been existed !!!  "                                                              )
    print( " If you want to overwrite the existing file directory, please input 'continue' in the terminal !! "     ) 
    print( " If you want to retain the existing file directory, please input 'stop' in the terminal to stop the "   ) 
    print( " simulation. Then you can reset the output dictionary in the input script file AMSS_NCKU_Input.py !!! " )
    print(                                                                                                          )
    ## 设定输入，是否覆盖原有文件夹
    while True:
        try:
            inputvalue = input()
            ## 如果同意覆盖当前目录文件目录，则退出计算，同时去除当前文件目录
            if ( inputvalue == "continue" ):
                print(                                  )
                print( " 继续计算 "                      )
                print( " Continue the calculation !!! " )
                print(                                  )
                break  # 输入合法，跳出循环
            ## 如果不同意覆盖当前目录文件目录，则退出计算，同时保留当前文件目录 
            elif ( inputvalue == "stop" ):
                print(                                 )
                print( " 退出计算 "                     )
                print( " Stop the calculation !!! "    )
                print(                                 )
                ## 退出整个程序
                sys.exit() 
            ## 如果用户没有按照要求输入，则要求用户重新输入 
            else:
                print(                                                    )
                print( " 请输入你的选择: 'continue' 或 'stop' !! "           )
                print( " Please input your choice !!! "                   )
                print( " Input 'continue' or 'stop' in the terminal !!! " )
        except ValueError:
            print(                                   )
            print( " 请输入你的选择: 'continue' 或 'stop' !! "           )
            print( " Please input your choice !!! "                   )
            print( " Input 'continue' or 'stop' in the terminal !!! " )
        
## 如果文件目录已经存在，则去除它
shutil.rmtree(File_directory, ignore_errors=True)

## 创建文件夹
os.mkdir(File_directory)

## 复制 python 输入文件到程序运行目录
shutil.copy("AMSS_NCKU_Input.py", File_directory)

# 生成文件目录，用来存放各类输出文件

output_directory = os.path.join(File_directory, "AMSS_NCKU_output")
os.mkdir(output_directory)

binary_results_directory = os.path.join(output_directory, input_data.Output_directory)
os.mkdir(binary_results_directory)

figure_directory = os.path.join(File_directory, "figure")
os.mkdir(figure_directory)

print(                                            )
print( " 生成文件目录完成 "                        )
print( " Output directionary has been generated " )
print(                                            )


##################################################################

## 输出相关的参数信息

import setup

## 输出相关的参数信息

setup.print_input_data( File_directory )
setup.generate_AMSSNCKU_input()

print(                                                                                           )
print( " 检查网格大小和分辨率是否满足要求 "                                                           )
print( " 如果网格大小和分辨率不合适，Ctrl+C 退出，调整网格层数和每层网格格点数目！！！ "                    )
print( " 如果网格大小和分辨率设定无误，按回车继续！！！ "                                               )
print(                                                                                           )
print( " Please check whether the grid boxes and their resolution is OK "                        )
print( " If the grid boxes and their resolution is not setting properly, press Ctrl+C to abort " )
print( " the simulation. Change the grid level structure and grid points !!! "                   )
print( " If the grid boxes and their resolution is appropriate, press Enter to continue !!!  "   )
inputvalue = input()  ## 设定一个输入（回车），以便程序下一步运行
print()

start_time = time.time()  # 记录开始时间

setup.print_puncture_information()


##################################################################

## 根据设定的参数生成 AMSS-NCKU 程序的输入文件

print(                                                                                   )
print( " 正在生成 AMSS-NCKU 程序的输入文件 "                                                )  
print( " Automatically generating the input parfile for AMSS-NCKU C++ program ABE.EXE. " ) 
print(                                                                                   ) 

## 根据格点信息生成 cgh 的相关输入

import numerical_grid        

numerical_grid.append_AMSSNCKU_cgh_input()  

print(                                                                                 )
print( " 完成 AMSS-NCKU 程序的输入文件 "                                                  )    
print( " 但是 AMSS-NCKU ABE 中 TwoPuncture 相关的输入要后几步运行后追加 "                   )
print( " The input parfile for AMSS-NCKU C++ executable file ABE has been generated. " )    
print( " However, the input relevant to TwoPuncture need to be appended later. "       )
print(                                                                                 )


##################################################################

## 画出初始的格点设置

print(                                                      )
print( " 正在画出初始格点的图像 "                              ) 
print( " Schematically plot the numerical grid structure. " ) 
print(                                                      )

numerical_grid.plot_initial_grid()


##################################################################

##  根据算法和参数生成 AMSS-NCKU 的宏文件

print(                                                                                   ) 
print( " 根据算法和参数生成 AMSS-NCKU 的宏文件 "                                             )       
print( " Automatically generating the macro file for AMSS-NCKU C++ executable file ABE " ) 
print( " (Based on the finite-difference numerical scheme) "                             )
print(                                                                                   )

import generate_macrodef

generate_macrodef.generate_macrodef_h()
print( " AMSS-NCKU 程序的宏文件 macrodef.h 已生成 "                       )   
print( " The macro file for AMSS-NCKU macrodef.h has been generated. " )
     
generate_macrodef.generate_macrodef_fh()
print( " AMSS-NCKU 程序的宏文件 macrodef.fh 已生成 "                       )  
print( " The macro file for AMSS-NCKU macrodef.fh has been generated. " )


##################################################################

##  根据用户的要求编译 AMSS-NCKU 程序

print(                                                         )
print( " 准备根据要求编译并运行 AMSS-NCKU 程序 "                   )
print( " Compiling the AMSS-NCKU code based on macro files "   )
print(                                                         )
#inputvalue = input()           
#print()

AMSS_NCKU_source_path = "AMSS_NCKU_source"
AMSS_NCKU_source_copy = os.path.join(File_directory, "AMSS_NCKU_source_copy")

###############################

## 如果 AMSS-NCKU 的源文件夹不存在，则创建它

# if not os.path.exists(destination_folder):
#     os.makedirs(destination_folder)

if not os.path.exists(AMSS_NCKU_source_path):
    os.makedirs(AMSS_NCKU_source_path)
    print( " 缺少 AMSS-NCKU 源文件，请将代码文件复制到 AMSS_NCKU_source 文件夹，按回车继续！！！ " )
    print( " The AMSS-NCKU source files are incomplete, please copy all source code files to the dictionary ./AMSS_NCKU_source. " )
    print( " Press Enter to continue!!! "                                                                                         )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input()
    
###############################

# 拷贝 AMSS-NCKU 的源文件，准备编译
shutil.copytree(AMSS_NCKU_source_path, AMSS_NCKU_source_copy)    

# 将整个src文件夹复制到dst位置
# shutil.copytree(src, dst)

# 将生成的宏文件拷贝进 AMSS-NCKU 的代码文件夹

macrodef_h_path  = os.path.join(File_directory, "macrodef.h") 
macrodef_fh_path = os.path.join(File_directory, "macrodef.fh") 

shutil.copy2(macrodef_h_path,  AMSS_NCKU_source_copy)
shutil.copy2(macrodef_fh_path, AMSS_NCKU_source_copy)

# 复制文件到目标位置
# shutil.copy2(source_file_path, target_location)
# shutil.copy2保留文件的元数据，如修改时间。如果你只想复制文件内容，不保留元数据，可以使用shutil.copy。

###############################

# 编译相关程序

import makefile_and_run

## 先切换到目标文件夹
os.chdir(AMSS_NCKU_source_copy)
 
## 编译 AMSS-NCKU 的主程序 ABE/ABEGPU
makefile_and_run.makefile_ABE()

## 如果初值方法选为 TwoPuncture，编译可执行文件 TwoPunctureABE 

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ): 
    makefile_and_run.makefile_TwoPunctureABE()
    
###########################

## 改变当前工作目录到上一级目录
os.chdir('..')
os.chdir('..')

print()

##################################################################

## 将 AMSS-NCKU 程序可执行文件 ABE 复制到程序运行的文件夹

if (input_data.GPU_Calculation == "no"):
    ABE_file = os.path.join(AMSS_NCKU_source_copy, "ABE")
elif (input_data.GPU_Calculation == "yes"):
    ABE_file = os.path.join(AMSS_NCKU_source_copy, "ABEGPU")

if not os.path.exists( ABE_file ):
    print(                                                                                                  )
    print( " 缺少 AMSS-NCKU 可执行文件 ABE/ABEGPU，请将 AMSS_NCKU_source 编译，编译完成后按回车继续！！！ "        )
    print( " Lack of AMSS-NCKU executable file ABE/ABEGPU, please recompile the ABE/ABEGPU in dictionary " )
    print( " AMSS_NCKU_source manually!!!   When recompile is finished, Press Enter to continue!!! "       )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input() 

## 复制可执行文件 ABE 到程序运行目录
shutil.copy2(ABE_file, output_directory)

###########################

## 如果初值方法选为 TwoPuncture，将 AMSS-NCKU 程序可执行文件 TwoPunctureABE 复制到程序运行的文件夹

TwoPuncture_file = os.path.join(AMSS_NCKU_source_copy, "TwoPunctureABE")

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):

    if not os.path.exists( TwoPuncture_file ):
        print(                                                                                                                  )
        print( " 缺少 AMSS-NCKU 可执行文件 TwoPunctureABE，请将 AMSS_NCKU_source 中的 TwoPunctureABE 编译，编译完成后按回车继续！！！ " )
        print( " Lack of AMSS-NCKU executable file TwoPunctureABE, please recompile the TwoPunctureABE in dictionary  "        ) 
        print( " AMSS_NCKU_source manually!!!   When recompile is finished, Press Enter to continue!!! "                       )
        inputvalue = input() 

    ## 复制可执行文件 TwoPunctureABE 到程序运行目录
    shutil.copy2(TwoPuncture_file, output_directory)

###########################


##################################################################

## 如果初值方法选为 TwoPuncture，生成 TwoPuncture 的输入文件

if (input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):

    print(                                                  )
    print( " 初始值类型选取为  Ansorg-Twopuncture"            )
    print( " Initial data is chosen as Ansorg-Twopuncture" )
    print(                                                 )
    
    print(                                                                                                         )
    print( " 正在生成 AMSS-NCKU TwoPuncture 程序的输入文件 "                                                          )  
    print( " Automatically generating the input parfile for AMSS-NCKU TwoPuncture executable file TwoPunctureABE " )
    print(                                                                                                         ) 
    
    import generate_TwoPuncture_input
    
    generate_TwoPuncture_input.generate_AMSSNCKU_TwoPuncture_input()
    
    print(                                                                                              )
    print( " 完成 AMSS-NCKU TwoPuncture 程序的输入文件 "                                                   )  
    print( " The input parfile for AMSS-NCKU TwoPuncture executable file TwoPunctureABE is generated. " )
    print(                                                                                              )
    
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
    
    makefile_and_run.run_TwoPunctureABE()
    
    ###########################
    
    ## 改变当前工作目录到上两级目录
    os.chdir('..')
    os.chdir('..')
    
##################################################################
    
## 根据 TwoPuncture 的运行结果来更新 puncture data

import renew_puncture_parameter
    
renew_puncture_parameter.append_AMSSNCKU_BSSN_input(File_directory, output_directory)


## 生成的 AMSS-NCKU 输入文件名
AMSS_NCKU_inputfile      = 'AMSS-NCKU.input'
AMSS_NCKU_inputfile_path = os.path.join(File_directory, AMSS_NCKU_inputfile)
 
## 复制并重命名文件
shutil.copy2( AMSS_NCKU_inputfile_path, os.path.join(output_directory, 'input.par') )


print(                                                                          )
print( " 成功将生成的 AMSS-NCKU 程序输入文件复制到运行目录 "                         )    
print( " Successfully copy all AMSS-NCKU input parfile to target dictionary. " )  
print(                                                                         )
    

##################################################################

##  启动 AMSS-NCKU 程序

print()
## print(" 准备启动 AMSS-NCKU 程序，按回车继续！！！ ")
## inputvalue = input()           
print()

## 先切换到目标文件夹
os.chdir( output_directory )

makefile_and_run.run_ABE()

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

print(                                                                          )
print( " 准备对 AMSS-NCKU 程序运行结果进行画图 "                                    )  
print( " Plotting the txt and binary results data in the AMSS-NCKU simulation " ) 
print(                                                                          )


import plot_xiaoqu
import plot_GW_strain_amplitude_xiaoqu

## 画出黑洞轨迹图
plot_xiaoqu.generate_puncture_orbit_plot(   binary_results_directory, figure_directory )
plot_xiaoqu.generate_puncture_orbit_plot3D( binary_results_directory, figure_directory )

## 画出黑洞间距随时间的变化图
plot_xiaoqu.generate_puncture_distence_plot( binary_results_directory, figure_directory )

## 画出引力波波形图
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_gravitational_wave_psi4_plot( binary_results_directory, figure_directory, i )
    plot_GW_strain_amplitude_xiaoqu.generate_gravitational_wave_amplitude_plot( binary_results_directory, figure_directory, i )

## 画出时空 ADM 质量的变化图
for i in range(input_data.Detector_Number):
    plot_xiaoqu.generate_ADMmass_plot( binary_results_directory, figure_directory, i )

## 画出哈密顿约束违反性况的变化图
for i in range(input_data.grid_level):
    plot_xiaoqu.generate_constraint_check_plot( binary_results_directory, figure_directory, i )

## 对储存的二进制数据画出图像
plot_xiaoqu.generate_binary_data_plot( binary_results_directory, figure_directory )

end_time = time.time()                # 记录结束时间
elapsed_time = end_time - start_time  # 计算运行时间

print(                                                 )
print( f" 程序运行时间 = {elapsed_time} 秒 "             )
print( f" This Program Cost = {elapsed_time} Seconds " )
print(                                                 )


##################################################################

print(                                                                                    )
print( " 程序顺利结束，谢谢您的使用 "                                                         )
print( " The AMSS-NCKU-Python simulation is successfully finished, thanks for using !!! " )
print(                                                                                    )

##################################################################


