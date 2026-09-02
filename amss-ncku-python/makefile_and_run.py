
##################################################################
##
## 这个文件设定了 AMSS-NCKU 编译和运行的相关命令
## 小曲
## 2025/01/24
##
##################################################################


import AMSS_NCKU_Input as input_data
import subprocess
import os


##################################################################



##################################################################

## 这个函数编译 AMSS-NCKU 主程序 ABE

def makefile_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                      )
        print( " 正在编译 AMSS-NCKU 程序 ABE/ABEGPU " ) 
        print(                                      )
    elif ( input_language == "English" ):  
        print(                                                        )
        print( " Compiling the AMSS-NCKU executable file ABE/ABEGPU " ) 
        print(                                                        )

    ## 编译命令
    if (input_data.GPU_Calculation == "no"):
        makefile_command  = "make -j4" + " ABE"
    elif (input_data.GPU_Calculation == "yes"):
        makefile_command  = "make -j4" + " ABEGPU"
    else:
        if ( input_language == "Chinese" ):
            print( " CPU/GPU 计算设置出错 " )
            print(                        )
        elif ( input_language == "English" ):
            print( " CPU/GPU numerical calculation setting is wrong " )
            print(                                                    )
 
    ## 使用subprocess.Popen来执行命令，并实时打印输出
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
 
    ## 循环读取输出并打印
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出
 
    ## 等待进程结束
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    if ( input_language == "Chinese" ):
        print(                              )
        print( " AMSS-NCKU 程序 ABE 编译完成" ) 
        print(                              )
    elif ( input_language == "English" ):
        print(                                                                  )
        print( " Compilation of the AMSS-NCKU executable file ABE is finished " ) 
        print(                                                                  )
    
    return
        
##################################################################



##################################################################

## 这个函数编译 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE

def makefile_TwoPunctureABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                          )
        print( " 正在编译 AMSS-NCKU 程序 TwoPunctureABE " ) 
        print(                                          )
    elif ( input_language == "English" ):
        print(                                                            )
        print( " Compiling the AMSS-NCKU executable file TwoPunctureABE " )
        print(                                                            )
    
    ## 编译命令
    makefile_command = "make" + " TwoPunctureABE"

    ## 使用subprocess.Popen来执行命令，并实时打印输出
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## 循环读取输出并打印
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出 
        
    ## 等待进程结束
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    if ( input_language == "Chinese" ):
        print(                                         )
        print( " AMSS-NCKU 程序 TwoPunctureABE 编译完成" ) 
        print(                                         )
    elif ( input_language == "English" ):   
        print(                                                                             ) 
        print( " Compilation of the AMSS-NCKU executable file TwoPunctureABE is finished " )
        print(                                                                             )
    
    return
    
##################################################################




##################################################################

## 这个函数编译 AMSS-NCKU 的后牛顿程序 PNOrbit

def makefile_PNOrbit( input_language ):

    if ( input_language == "Chinese" ):
        print(                                    )
        print( " 正在编译 AMSS-NCKU 程序 PNOrbit " ) 
        print(                                    )
    elif ( input_language == "English" ):
        print(                                                      )
        print( " Compiling the AMSS-NCKU executable file PNOrbit " )
        print(                                                      )
    
    ## 编译命令
    makefile_command = "make" + " PNOrbit"

    ## 使用subprocess.Popen来执行命令，并实时打印输出
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## 循环读取输出并打印
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出 
        
    ## 等待进程结束
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    if ( input_language == "Chinese" ):
        print(                                   )
        print( " AMSS-NCKU 程序 PNOrbit 编译完成" ) 
        print(                                   )
    elif ( input_language == "English" ):   
        print(                                                                       ) 
        print( " Compilation of the AMSS-NCKU executable file PNOrbit is finished " )
        print(                                                                       )
    
    return
    
##################################################################




##################################################################

## 这个函数运行 AMSS-NCKU 主程序 ABE

def run_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                        )
        print( " 正在运行 AMSS-NCKU 主程序 ABE/ABEGPU " ) 
        print(                                        )
    elif ( input_language == "English" ): 
        print(                                                      )    
        print( " Running the AMSS-NCKU executable file ABE/ABEGPU " ) 
        print(                                                      )

    ## 定义要运行的命令，要使用 str 将其它转换为字符串
    
    if (input_data.GPU_Calculation == "no"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABE"
        mpi_command_outfile = "ABE_out.log"
    elif (input_data.GPU_Calculation == "yes"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABEGPU"
        mpi_command_outfile = "ABEGPU_out.log"
 
    ## 使用subprocess.Popen来执行命令，并实时打印输出
    mpi_process = subprocess.Popen(mpi_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
 
    ## 将 ABE 的运行结果写入文件中
    with open(mpi_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印
        for line in mpi_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # 将行写入文件
            file0.flush()        # 确保每行都被立即写入文件，可选            
    file0.close()
 
    ## 等待进程结束
    mpi_return_code = mpi_process.wait()
    
    if ( input_language == "Chinese" ):
        print(                                        )
        print( " AMSS-NCKU 主程序 ABE/ABEGPU 运行结束 " ) 
        print(                                        )
    elif ( input_language == "English" ): 
        print(                                           )
        print( " The ABE/ABEGPU simulation is finished " ) 
        print(                                           )
    
    return

##################################################################



##################################################################

## 这个函数运行 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE

def run_TwoPunctureABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                          )
        print( " 正在运行 AMSS-NCKU 程序 TwoPunctureABE " )
        print(                                          ) 
    elif ( input_language == "English" ): 
        print(                                                          ) 
        print( " Running the AMSS-NCKU executable file TwoPunctureABE " ) 
        print(                                                          )
    
    ## 定义要运行的命令，要使用 str 将其它转换为字符串
    
    TwoPuncture_command         = "./TwoPunctureABE"
    TwoPuncture_command_outfile = "TwoPunctureABE_out.log"

    ## 设置 OpenMP 多线程数目，线程数取输入文件 AMSS_NCKU_Input.py 中设定的 OMP_processes
    ## 通过环境变量 OMP_NUM_THREADS 传给 TwoPunctureABE，无需修改 C 程序
    TwoPuncture_env = os.environ.copy()
    TwoPuncture_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## 使用subprocess.Popen来执行命令，并实时打印输出
    
    TwoPuncture_process = subprocess.Popen(TwoPuncture_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=TwoPuncture_env)
 
    ## 将 TwoPunctureABE 的运行结果写入文件中
    
    with open(TwoPuncture_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印
        for line in TwoPuncture_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # 将行写入文件
            file0.flush()        # 确保每行都被立即写入文件，可选                 
    file0.close()
 
    ## 等待进程结束
    TwoPuncture_command_return_code = TwoPuncture_process.wait()
    
    if ( input_language == "Chinese" ):
        print(                                          )
        print( " AMSS-NCKU 程序 TwoPunctureABE 运行结束 " )
        print(                                          ) 
    elif ( input_language == "English" ):
        print(                                               )
        print( " The TwoPunctureABE simulation is finished " ) 
        print(                                               )
    
    return

##################################################################



##################################################################

## 这个函数运行 AMSS-NCKU 的后牛顿程序 PNOrbit

def run_PNOrbit( input_language ):

    if ( input_language == "Chinese" ):
        print(                                   )
        print( " 正在运行 AMSS-NCKU 程序 PNOrbit " )
        print(                                   ) 
    elif ( input_language == "English" ): 
        print(                                                   ) 
        print( " Running the AMSS-NCKU executable file PNOrbit " ) 
        print(                                                   )
    
    ## 定义要运行的命令，要使用 str 将其它转换为字符串
    
    PNOrbit_command         = "./PNOrbit"
    PNOrbit_command_outfile = "PNOrbit_out.log"

    ## 设置 OpenMP 多线程数目，线程数取输入文件 AMSS_NCKU_Input.py 中设定的 OMP_processes
    ## 通过环境变量 OMP_NUM_THREADS 传给 PNOrbit，无需修改 C 程序
    PNOrbit_env = os.environ.copy()
    PNOrbit_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## 使用subprocess.Popen来执行命令，并实时打印输出
    
    PNOrbit_process = subprocess.Popen(PNOrbit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=PNOrbit_env)
 
    ## 将 PNOrbit 的运行结果写入文件中
    
    with open(PNOrbit_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印
        for line in PNOrbit_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # 将行写入文件
            file0.flush()        # 确保每行都被立即写入文件，可选                 
    file0.close()
 
    ## 等待进程结束
    PNOrbit_command_return_code = PNOrbit_process.wait()
    
    if ( input_language == "Chinese" ):
        print(                                   )
        print( " AMSS-NCKU 程序 PNOrbit 运行结束 " )
        print(                                   ) 
    elif ( input_language == "English" ):
        print(                                        )
        print( " The PNOrbit simulation is finished " ) 
        print(                                        )
    
    return

##################################################################
    
