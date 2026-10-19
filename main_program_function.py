
##################################################################
##
## 这个文件设定了 AMSS-NCKU Python 主程序的一些特殊功能 | This file sets up some special functions of the AMSS-NCKU Python main program
## 小曲 | Xiaoqu
## 2026/02/15
##
##################################################################

import sys
import os
import generate_macrodef
import AMSS_NCKU_Input as input_data

##################################################################

## 选择终端的输出语言 | Choose the output language of the terminal

def setting_language():

    print(                                                      )
    print( " Please Choosing the Language for terminal output " )
    print( " Input 'English' or 'Chinese' "                     )
    print(                                                      )

    ## 如果属于语言不正确，则重新输入，直到用户输入正确的语言 | If the input language is incorrect, re-enter it until the user inputs a valid language
    while True:
        try: 
            ## 输入选择的语言 | Input the chosen language
            input_language = input() 
        
            ## 如果选择汉语，可以跳出循环 | If Chinese is chosen, break out of the loop
            if ( input_language == "Chinese" ):
                print(                                       )
                print( " 选择的屏幕输出语言为 ", input_language )         
                print(                                       )
                break  # 输入合法，跳出循环 | Valid input, break out of the loop
            ## 如果选择英语，可以跳出循环 | If English is chosen, break out of the loop
            elif ( input_language == "English" ):
                print(                                        )
                print( " Chosen Language is ", input_language )         
                print(                                        )
                break  # 输入合法，跳出循环 | Valid input, break out of the loop
            ## 如果用户没有按照要求输入，则要求用户重新输入 | If the user does not input as required, ask the user to re-enter
            else:
                print(                                                              )
                print( " Please Choosing the correct Language for terminal output " )
                print( " Input 'English' or 'Chinese' "                             )
                print( " Other languages are not supproted "                        )
                print(                                                              )
        except ValueError:
            print(                                                              )
            print( " Please Choosing the correct Language for terminal output " )
            print( " Input 'English' or 'Chinese' "                             )
            print( " Other languages are not supproted "                        )
            print(                                                              )
            
    return input_language

##################################################################



##################################################################

## 当设定的文件夹存在时 | When the configured folder exists
## 这个函数用来设定是否覆盖上一次计算 | This function sets whether to overwrite the previous calculation

def continue_or_stop( File_directory, input_language ):

    ## 如果设定的输出目录存在，则根据用户的选择看是否继续计算 | If the configured output directory exists, decide whether to continue the calculation according to the user's choice
    if os.path.exists(File_directory):
    
        if ( input_language == "Chinese" ):
            print(                                                                                         )
            print( " 设定的输出目录存在，是否同意覆盖当前目录？ "                                                 )
            print( " 如果同意覆当前目录，请输入'continue'继续计算 "                                              )
            print( " 如果不同意覆当前目录，请输入'stop'退出计算，并在输入文件 AMSS_NCKU_Input.py 中重新设置输出目录 "  )
            print(                                                                                         )
        elif ( input_language == "English" ):    
            print(                                                                                                          )
            print( " Output dictionary has been existed !!!  "                                                              )
            print(                                                                                                          )
            print( " If you want to overwrite the existing file directory, please input 'continue' in the terminal !! "     ) 
            print(                                                                                                          )
            print( " If you want to retain the existing file directory, please input 'stop' in the terminal to stop the "   ) 
            print( " simulation. Then you can reset the output dictionary in the input script file AMSS_NCKU_Input.py !!! " )
            print(                                                                                                          )
        
        ## 设定输入，是否覆盖原有文件夹 | Set the input: whether to overwrite the existing folder
        
        while True:
            try:
                inputvalue = input()
                ## 如果同意覆盖当前目录文件目录，则退出计算，同时去除当前文件目录 | If overwriting the current file directory is agreed, exit the calculation and remove the current file directory
                if ( inputvalue == "continue" ):
                    if ( input_language == "Chinese" ):
                        print(             )
                        print( " 继续计算 " )
                        print(             )
                    elif ( input_language == "English" ):
                        print(                                  )
                        print( " Continue the calculation !!! " )
                        print(                                  )
                    break  # 输入合法，跳出循环 | Valid input, break out of the loop
                ## 如果不同意覆盖当前目录文件目录，则退出计算，同时保留当前文件目录 | If overwriting the current file directory is not agreed, exit the calculation and keep the current file directory 
                elif ( inputvalue == "stop" ):
                    if ( input_language == "Chinese" ):
                        print(             )
                        print( " 退出计算 " )
                        print(             )
                    elif ( input_language == "English" ):
                        print(                                 )
                        print( " Stop the calculation !!! "    )
                        print(                                 )
                    ## 退出整个程序 | Exit the whole program
                    sys.exit() 
                ## 如果用户没有按照要求输入，则要求用户重新输入 | If the user does not input as required, ask the user to re-enter 
                else:
                    if ( input_language == "Chinese" ):
                        print(                                          )
                        print( " 请输入你的选择: 'continue' 或 'stop' !! " )
                        print(                                          )
                    elif ( input_language == "English" ):
                        print(                                                    )
                        print( " Please input your choice !!! "                   )
                        print( " Input 'continue' or 'stop' in the terminal !!! " )
                        print(                                                    )
            except ValueError:
                if ( input_language == "Chinese" ):
                    print(                                            )
                    print( " 请输入你的选择: 'continue' 或 'stop' !! " )
                    print(                                            )
                elif ( input_language == "English" ):
                    print(                                                    )
                    print( " Please input your choice !!! "                   )
                    print( " Input 'continue' or 'stop' in the terminal !!! " )
                    print(                                                    )
                    
    return


##################################################################



##################################################################

##  这个函数根据算法和参数生成 AMSS-NCKU 的宏文件 | This function generates the AMSS-NCKU macro file from the algorithm and parameters

def generate_macrodef_file( input_language ):

    if ( input_language == "Chinese" ):
        print(                                       ) 
        print( " 根据算法和参数生成 AMSS-NCKU 的宏文件 " ) 
        print(                                       ) 
    elif ( input_language == "English" ): 
        print(                                                                               )     
        print( " Automatically generating the macro file for AMSS-NCKU C++ executable file " ) 
        print( " (Based on the finite-difference numerical scheme) "                         )
        print(                                                                               )

    generate_macrodef.generate_macrodef_h()
    
    if ( input_language == "Chinese" ):
        print( " AMSS-NCKU 程序的宏文件 macrodef.h 已生成 " )  
    elif ( input_language == "English" ):  
        print( " The macro file for AMSS-NCKU macrodef.h has been generated." )
     
    generate_macrodef.generate_macrodef_fh()
    
    if ( input_language == "Chinese" ):
        print( " AMSS-NCKU 程序的宏文件 macrodef.fh 已生成 " )  
    elif ( input_language == "English" ):
        print( " The macro file for AMSS-NCKU macrodef.fh has been generated." )
        
    return
    
##################################################################



##################################################################

##  这个函数根据输入文件中选择的计算方式（CPU 或 GPU） | This function determines, according to the calculation mode (CPU or GPU) chosen in the input file,
##  以及 GPU 内层代码的选择（flash 或 pro） | and the choice of the GPU inner code (flash or pro),
##  确定编译时所使用的源代码文件夹 | the source-code folder used for compilation

def select_source_path( input_language ):

    AMSS_NCKU_shared_path = "AMSS_NCKU_source"

    if ( input_data.GPU_Calculation == "no" ):
        ## 选择 CPU 计算时使用 AMSS_NCKU_source/ABE_CPU | Use AMSS_NCKU_source/ABE_CPU when CPU calculation is chosen
        AMSS_NCKU_source_path = os.path.join( AMSS_NCKU_shared_path, "ABE_CPU" )
    elif ( input_data.GPU_Calculation == "yes" ):
        ## 选择 GPU 计算时根据 GPU_Code_Choice 选择源代码文件夹 | When GPU calculation is chosen, select the source-code folder according to GPU_Code_Choice
        if ( input_data.GPU_Code_Choice == "flash" ):
            AMSS_NCKU_source_path = os.path.join( AMSS_NCKU_shared_path, "ABE_GPU_flash" )
        elif ( input_data.GPU_Code_Choice == "pro" ):
            AMSS_NCKU_source_path = os.path.join( AMSS_NCKU_shared_path, "ABE_GPU_pro" )
        else:
            if ( input_language == "Chinese" ):
                print( " GPU 代码选择设置出错，请将 AMSS_NCKU_Input.py 中的 GPU_Code_Choice 设为 flash 或 pro " )
            elif ( input_language == "English" ):
                print( " GPU code choice setting is wrong, please set GPU_Code_Choice to flash or pro in AMSS_NCKU_Input.py " )
            raise SystemExit(1)
    else:
        if ( input_language == "Chinese" ):
            print( " CPU/GPU 计算设置出错，请将 AMSS_NCKU_Input.py 中的 GPU_Calculation 设为 yes 或 no " )
        elif ( input_language == "English" ):
            print( " CPU/GPU numerical calculation setting is wrong, please set GPU_Calculation to yes or no " )
        raise SystemExit(1)

    return AMSS_NCKU_source_path
    
##################################################################

