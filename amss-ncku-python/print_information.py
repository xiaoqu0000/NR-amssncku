
##################################################################
##
## 该文件定义程序的介绍
## 小曲
## 2025/02/07
##
##################################################################



##################################################################

## 这个函数用来输出整个程序的介绍

def print_program_introduction_Chinese():

    print(                                                                                                  )
    print( "------------------------------------------------------------------------------------------"     )  
    print(                                                                                                  )
    print( "     数值相对论计算程序 AMSS-NCKU  "                                                              )
    print(                                                                                                 )
    print( "     原程序作者： 曹周键等 "                                                                      )
    print( "     程序 Python 接口作者：小曲 "                                                                 )
    print(                                                                                                 )
    print( "     AMSS-NCKU 是一个数值相对论程序  "                                                            )
    print( "     本程序用来对双黑洞合并过程进行数值模拟，通过数值求解爱因斯坦方程得到双黑洞合并过程中引力场随时间的演化，"  )
    print( " 从而得到黑洞的轨迹和释放引力波的信息。"                                                             )
    print(                                                                                                 )
    print( "     在数值方法上，本程序用使用有限差分方法对双黑洞合并过程进行数值模拟。程序中可以选择的差分方法有：2 阶 "    )
    print( " 差分、4 阶差分、6 阶差分、8 阶差分 "                                                              )
    print( "     程序中可以选择的微分方程为：BSSN 方程、Z4C 方程、BSSN 方程耦合标量场、BSSN 方程耦合电磁场 "          ) 
    print( "     程序中可以选择的网格类型为：方形网格、最外层带球壳的 shell patch 网格 "                            )
    print(                                                                                                 )
    print( "     除此之外，本程序还实现了 CPU 和 GPU 的混合运算"                                                 )
    print(                                                                                                 )
    print( "------------------------------------------------------------------------------------------"    ) 
    print(                                                                                                 )
    
    return
    
##################################################################    

## 这个函数用来输出整个程序的介绍
## 这是同一份函数，用英语输出 
    
def print_program_introduction_English():

    print(                                                                                                          )
    print( "------------------------------------------------------------------------------------------"             )  
    print(                                                                                                          )
    print( "     Numerical Relativity AMSS-NCKU  "                                                                  )
    print(                                                                                                          )
    print( "     Author of AMSS-NCKU Code: Zhou-Jian Cao et al. "                                                   )
    print( "     Author of AMSS-NCKU Python Interface: Xiao Qu "                                                    )
    print(                                                                                                          )
    print( "     AMSS-NCKU is an open source numerical relativity code "                                            )
    print( "     It can be used to simulate the dynamical evolution on mergering process of black hole systems, "   )
    print( " calculating the variation of gravitational field, black holes' trajectories, and gravitational wave "  )
    print( " emissions through directly solving the Einstein field equations  "                                     )
    print(                                                                                                          )
    print( "     This AMSS-NCKU code uses the finite-difference method to evaluate the numerical simulation. The "  )
    print( " finite-difference schemes can be chosen as: 2nd order, 4th order, 6th order, 8th order. "              )
    print( "     The computation equation form in AMSS-NCKU code can be chosen as: BSSN equations, Z4C equations, " )
    print( " BSSN equations coupled with scalars (in f(R) theory), BSSN equations coupled with electromagnetic "    )
    print( " fields. "                                                                                              ) 
    print( "     The numerical grid system in this code includes: patch AMR grid, shell-patch AMR grid. "           )
    print(                                                                                                          )
    print( "     Furthermore, This code has fulfilled the CPU and GPU hybrid calculation. "                         )
    print(                                                                                                          )
    print( "------------------------------------------------------------------------------------------"             ) 
    print(                                                                                                          )
    
    return
    
##################################################################



##################################################################

## 该函数用于在屏幕上输出程序开始的提示

def print_begin_program( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                                  )
        print( " 计算即将开始，请确认在 AMSS_NCKU_Input.py 中设置了正确的参数，按回车继续！！！  " )
        print( " 如果输入参数没有设置好，Ctrl+C 退出，调整 AMSS_NCKU_Input.py 中的输入参数！！！ " )
        print(                                                                                  )
    elif ( input_language == "English" ):
        print(                                                                                                       )
        print( " Simulation will be started, please confirm you have set the correct parameters in the script file " )
        print( " AMSS_NCKU_Input.py "                                                                                )
        print(                                                                                                       )
        print( " If parameters have been set correctly, press Enter to continue !!!  "                               )
        print(                                                                                                       )
        print( " If you have not set parameters, press Ctrl+C to abort the simulation and adjust the parameters "    )
        print( " in script file AMSS_NCKU_Input.py !!! "                                                             )
        print(                                                                                                       )
    
    return 
    
##################################################################  


################################################################## 

## 该函数用于在屏幕上输出文件目录已生成

def print_make_directory( input_language ):

    if ( input_language == "Chinese" ):
        print(                    )
        print( " 生成文件目录完成 "  )
        print(                    )
    elif ( input_language == "English" ):
        print(                                         )
        print( " Output directory has been generated " )
        print(                                         ) 
    
    return
    
################################################################## 



##################################################################  

## 该函数用于在在屏幕上输出，并指引用户判断是否需要开始运算

def print_whether_grid_is_satisfied( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                       )
        print( " 检查网格大小和分辨率是否满足要求 "                                       )
        print( " 如果网格大小和分辨率不合适，Ctrl+C 退出，调整网格层数和每层网格格点数目！！！ " )
        print( " 如果网格大小和分辨率设定无误，按回车继续！！！ "                            )
        print(                                                                       )
    elif ( input_language == "English" ):
        print(                                                                                             )
        print( " Please check whether the grid boxes and their resolution are OK "                         )
        print(                                                                                             )
        print( " If the grid boxes and their resolutions are not setting properly, press Ctrl+C to abort " )
        print( " the simulation. Change the grid level structure and grid points !!! "                     )
        print(                                                                                             )
        print( " If the grid boxes and their resolutions are appropriate, press Enter to continue !!!  "   )
        print(                                                                                             )
    
    return

##################################################################  




################################################################## 

## 该函数用于输出开始编译的消息

def print_compile_AMSS_NCKU( input_language ): 

    if ( input_language == "Chinese" ):
        print(                                       )
        print( " 准备根据要求编译并运行 AMSS-NCKU 程序 " )
        print(                                       )
    elif ( input_language == "English" ):
        print(                                                       )
        print( " Compiling the AMSS-NCKU code based on macro files " )
        print(                                                       )
    #inputvalue = input()           
    #print()

    return

################################################################## 



################################################################## 

## 该函数用于输出编译 AMSS_NCKU 出错的信息

def print_compile_AMSS_NCKU_debug( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                                )
        print( " 缺少 AMSS-NCKU 源文件，请将代码文件复制到 AMSS_NCKU_source 文件夹，按回车继续！！！ " )
        print(                                                                                )
    elif ( input_language == "English" ):
        print(                                                                             )
        print( " The AMSS-NCKU source files are incomplete !! "                            )
        print( " Please copy all source code files to the dictionary ./AMSS_NCKU_source. " )
        print( " Press Enter to continue!!! "                                              )
        print(                                                                             )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input()



##################################################################   

## 该函数用输出成功拷贝参数文件

def copy_input_parfile( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                 )
        print( " 成功将生成的 AMSS-NCKU 程序输入文件复制到运行目录 " )  
        print(                                                 )
    elif ( input_language == "English" ):  
        print(                                                                         )
        print( " Successfully copy all AMSS-NCKU input parfile to target dictionary. " )  
        print(                                                                         )
        
    return
    
##################################################################



##################################################################   

## 该函数用于在 AMSS-NCKU 的 ABE/ABEGPU 开始启动时输出文字信息

def print_begin_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                             )
        print( " 准备启动 AMSS-NCKU 的可执行程序 ABE/ABEGPU " )          
        print(                                             )
    elif ( input_language == "English" ):
        print(                                                                                    )
        print( " Get ready to launch the ABE/ABEGPU executable file in the AMSS-NCKU simulation " )          
        print(                                                                                    )
        
    return
    
##################################################################



##################################################################   

## 该函数用于输出 ABE/ABEGPU 不存在时的报错信息

def print_no_ABE_error( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                                         )
        print( " 缺少 AMSS-NCKU 可执行文件 ABE/ABEGPU，请将 AMSS_NCKU_source 编译，编译完成后按回车继续！！！ " )
        print(                                                                                         )
    elif ( input_language == "English" ):
        print(                                                                                                 )
        print( " Lack of AMSS-NCKU executable file ABE/ABEGPU, please recompile the ABE/ABEGPU in dictionary " )
        print( " AMSS_NCKU_source manually!!! "                                                                )
        print( " When recompile is finished, Press Enter to continue!!! "                                      )
        print(                                                                                                 )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input() 
        
    return
    
##################################################################




##################################################################

## 该函数用于输出生成 TwoPuncture 输入文件的消息

def print_generate_TwoPunture_input( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                      )
        print( " 初始值类型选取为  Ansorg-Twopuncture"                 )
        print(                                                      )
        print( " 正在生成 AMSS-NCKU 可执行文件 TwoPuncture 的输入文件 "  )
        print(                                                      )
    elif ( input_language == "English" ):
        print(                                                                                             ) 
        print( " Initial data is chosen as Ansorg-Twopuncture"                                             )
        print(                                                                                             ) 
        print( " Automatically generating the input parfile for AMSS-NCKU executable file TwoPunctureABE " )
        print(                                                                                             ) 
        
    return


##################################################################



##################################################################

## 该函数用于输出完成 TwoPuncture 输入文件的消息

def print_finish_TwoPunture_input( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                  )
        print( " 完成 AMSS-NCKU 可执行文件 TwoPuncture 的输入文件 " ) 
        print(                                                  )
    elif ( input_language == "English" ): 
        print(                                                                                  )
        print( " The input parfile for AMSS-NCKU executable file TwoPunctureABE is generated. " )
        print(                                                                                  )
        
    return

##################################################################



##################################################################

## 该函数用于输出 TwoPunture 不存在时的报错信息

def print_no_TwoPunture_error( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                        )
        print( " 缺少 AMSS-NCKU 可执行文件 TwoPunctureABE "                              )
        print( "请将 AMSS_NCKU_source 中的 TwoPunctureABE 编译，编译完成后按回车继续！！！ " )
        print(                                                                        )
    elif ( input_language == "English" ):
        print(                                                                                                         )
        print( " Lack of AMSS-NCKU executable file TwoPunctureABE, please recompile the TwoPunctureABE in dictionary " ) 
        print( " AMSS_NCKU_source manually!!!  "                                                                       )
        print( " When recompile is finished, Press Enter to continue!!! "                                              )
        print(                                                                                                         )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input()

    return
    
##################################################################




##################################################################

## 该函数用于输出生成 PNOrbit 输入文件的消息

def print_generate_PNOrbit_input( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                      )
        print( " 设置 puncture 数据的方式选取为  Automatically-BBH"    )
        print(                                                      )
        print( " 正在生成 AMSS-NCKU 可执行文件 PNOrbit 程序的输入文件 "  )
        print(                                                     )
    elif ( input_language == "English" ):
        print(                                                                                      ) 
        print( " Puncture data setteing is chosen as Automatically-BBH"                             )
        print(                                                                                      ) 
        print( " Automatically generating the input parfile for AMSS-NCKU executable file PNOrbit " )
        print(                                                                                      ) 
        
    return


##################################################################



##################################################################

## 该函数用于输出完成 PNOrbit 输入文件的消息

def print_finish_PNOrbit_input( input_language ):

    if ( input_language == "Chinese" ):
        print(                                              )
        print( " 完成 AMSS-NCKU 可执行文件 PNOrbit 的输入文件 " ) 
        print(                                              )
    elif ( input_language == "English" ): 
        print(                                                                           )
        print( " The input parfile for AMSS-NCKU executable file PNOrbit is generated. " )
        print(                                                                           )
        
    return

##################################################################



##################################################################

## 该函数用于输出 PNOrbit 不存在时的报错信息

def print_no_PNOrbit_error( input_language ):

    if ( input_language == "Chinese" ):
        print(                                                                  )
        print( " 缺少 AMSS-NCKU 可执行文件 PNOrbit "                              )
        print( "请将 AMSS_NCKU_source 中的 PNOrbit 编译，编译完成后按回车继续！！！ " )
        print(                                                                 )
    elif ( input_language == "English" ):
        print(                                                                                           )
        print( " Lack of AMSS-NCKU executable file PNOrbit, please recompile the PNOrbit in dictionary " ) 
        print( " AMSS_NCKU_source manually!!!  "                                                         )
        print( " When recompile is finished, Press Enter to continue!!! "                                )
        print(                                                                                           )
    ## 设定一个输入（回车），以便程序下一步运行
    inputvalue = input()

    return
    
##################################################################



    
##################################################################   

## 该函数用于在开始画图时输出文字信息

def print_begin_plot( input_language ):

    if ( input_language == "Chinese" ):
        print(                                       )
        print( " 准备对 AMSS-NCKU 程序运行结果进行画图 " )
        print(                                       )
    elif ( input_language == "English" ):  
        print(                                                                          )
        print( " Plotting the txt and binary results data in the AMSS-NCKU simulation " ) 
        print(                                                                          )

    return
    
################################################################## 



##################################################################

## 该函数用于输出计算所用的时间

def print_time_cost( elapsed_time, input_language ):

    if ( input_language == "Chinese" ):
        print(                                          )
        print( f" 程序运行共花费时间 = {elapsed_time} 秒 " )
        print(                                          )
    elif ( input_language == "English" ):
        print(                                                                         )
        print( f" The computer time cost in this simulation = {elapsed_time} Seconds " )
        print(                                                                         )

    return
    
##################################################################



##################################################################

## 该函数用于输出程序结束的提示

def print_program_end( input_language ):

    if ( input_language == "Chinese" ):
        print(                            )
        print( " 程序顺利结束，谢谢您的使用 " )
        print(                            )
    elif ( input_language == "English" ):
        print(                                                                                    )
        print( " The AMSS-NCKU-Python simulation is successfully finished, thanks for using !!! " )
        print(                                                                                    )
        
    return
    
##################################################################

