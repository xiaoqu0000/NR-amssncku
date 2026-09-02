
##################################################################
##
## 这个文件设定了 AMSS-NCKU Python 主程序的一些特殊功能
## 小曲
## 2026/02/15
##
##################################################################

import sys
import os
import generate_macrodef

##################################################################

## 选择终端的输出语言

def setting_language():

    print(                                                      )
    print( " Please Choosing the Language for terminal output " )
    print( " Input 'English' or 'Chinese' "                     )
    print(                                                      )

    ## 如果属于语言不正确，则重新输入，直到用户输入正确的语言
    while True:
        try: 
            ## 输入选择的语言
            input_language = input() 
        
            ## 如果选择汉语，可以跳出循环
            if ( input_language == "Chinese" ):
                print(                                       )
                print( " 选择的屏幕输出语言为 ", input_language )         
                print(                                       )
                break  # 输入合法，跳出循环
            ## 如果选择英语，可以跳出循环
            elif ( input_language == "English" ):
                print(                                        )
                print( " Chosen Language is ", input_language )         
                print(                                        )
                break  # 输入合法，跳出循环
            ## 如果用户没有按照要求输入，则要求用户重新输入
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

## 当设定的文件夹存在时
## 这个函数用来设定是否覆盖上一次计算

def continue_or_stop( File_directory, input_language ):

    ## 如果设定的输出目录存在，则根据用户的选择看是否继续计算
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
        
        ## 设定输入，是否覆盖原有文件夹
        
        while True:
            try:
                inputvalue = input()
                ## 如果同意覆盖当前目录文件目录，则退出计算，同时去除当前文件目录
                if ( inputvalue == "continue" ):
                    if ( input_language == "Chinese" ):
                        print(             )
                        print( " 继续计算 " )
                        print(             )
                    elif ( input_language == "English" ):
                        print(                                  )
                        print( " Continue the calculation !!! " )
                        print(                                  )
                    break  # 输入合法，跳出循环
                ## 如果不同意覆盖当前目录文件目录，则退出计算，同时保留当前文件目录 
                elif ( inputvalue == "stop" ):
                    if ( input_language == "Chinese" ):
                        print(             )
                        print( " 退出计算 " )
                        print(             )
                    elif ( input_language == "English" ):
                        print(                                 )
                        print( " Stop the calculation !!! "    )
                        print(                                 )
                    ## 退出整个程序
                    sys.exit() 
                ## 如果用户没有按照要求输入，则要求用户重新输入 
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

##  这个函数根据算法和参数生成 AMSS-NCKU 的宏文件

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

