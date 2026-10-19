
##################################################################
##
## 这个文件设定了 AMSS-NCKU 编译和运行的相关命令 | This file sets the commands related to compiling and running AMSS-NCKU
## 小曲 | Xiaoqu
## 2025/01/24
##
##################################################################


import AMSS_NCKU_Input as input_data
import subprocess
import os


##################################################################



##################################################################

## 这个函数编译 AMSS-NCKU 主程序 ABE | This function compiles the AMSS-NCKU main program ABE

def makefile_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                      )
        print( " 正在编译 AMSS-NCKU 程序 ABE/ABEGPU " ) 
        print(                                      )
    elif ( input_language == "English" ):  
        print(                                                        )
        print( " Compiling the AMSS-NCKU executable file ABE/ABEGPU " ) 
        print(                                                        )

    ## 编译命令 | Compilation command
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
 
    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
 
    ## 循环读取输出并打印 | Read the output in a loop and print it
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出 | Print the output in real time
 
    ## 等待进程结束 | Wait for the process to finish
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

## 这个函数编译 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE | This function compiles the AMSS-NCKU TwoPuncture program TwoPunctureABE

def makefile_TwoPunctureABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                          )
        print( " 正在编译 AMSS-NCKU 程序 TwoPunctureABE " ) 
        print(                                          )
    elif ( input_language == "English" ):
        print(                                                            )
        print( " Compiling the AMSS-NCKU executable file TwoPunctureABE " )
        print(                                                            )
    
    ## 编译命令 | Compilation command
    makefile_command = "make" + " TwoPunctureABE"

    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## 循环读取输出并打印 | Read the output in a loop and print it
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出 | Print the output in real time 
        
    ## 等待进程结束 | Wait for the process to finish
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

## 这个函数编译 AMSS-NCKU 的后牛顿程序 PNOrbit | This function compiles the AMSS-NCKU post-Newtonian program PNOrbit

def makefile_PNOrbit( input_language ):

    if ( input_language == "Chinese" ):
        print(                                    )
        print( " 正在编译 AMSS-NCKU 程序 PNOrbit " ) 
        print(                                    )
    elif ( input_language == "English" ):
        print(                                                      )
        print( " Compiling the AMSS-NCKU executable file PNOrbit " )
        print(                                                      )
    
    ## 编译命令 | Compilation command
    makefile_command = "make" + " PNOrbit"

    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## 循环读取输出并打印 | Read the output in a loop and print it
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出 | Print the output in real time 
        
    ## 等待进程结束 | Wait for the process to finish
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

## 这个函数运行 AMSS-NCKU 主程序 ABE | This function runs the AMSS-NCKU main program ABE

def run_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                        )
        print( " 正在运行 AMSS-NCKU 主程序 ABE/ABEGPU " ) 
        print(                                        )
    elif ( input_language == "English" ): 
        print(                                                      )    
        print( " Running the AMSS-NCKU executable file ABE/ABEGPU " ) 
        print(                                                      )

    ## 定义要运行的命令，要使用 str 将其它转换为字符串 | Define the command to run; use str to convert it into a string
    
    if (input_data.GPU_Calculation == "no"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABE"
        mpi_command_outfile = "ABE_out.log"
    elif (input_data.GPU_Calculation == "yes"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABEGPU"
        mpi_command_outfile = "ABEGPU_out.log"
 
    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    mpi_process = subprocess.Popen(mpi_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
 
    ## 将 ABE 的运行结果写入文件中 | Write the running results of ABE to a file
    with open(mpi_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印 | Read the output in a loop and print it
        for line in mpi_process.stdout:
            print(line, end='')  # 实时打印输出 | Print the output in real time
            file0.write(line)    # 将行写入文件 | Write the line to the file
            file0.flush()        # 确保每行都被立即写入文件，可选 | Ensure that each line is written to the file immediately (optional)            
    file0.close()
 
    ## 等待进程结束 | Wait for the process to finish
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

## 这个函数运行 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE | This function runs the AMSS-NCKU TwoPuncture program TwoPunctureABE

def run_TwoPunctureABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                          )
        print( " 正在运行 AMSS-NCKU 程序 TwoPunctureABE " )
        print(                                          ) 
    elif ( input_language == "English" ): 
        print(                                                          ) 
        print( " Running the AMSS-NCKU executable file TwoPunctureABE " ) 
        print(                                                          )
    
    ## 定义要运行的命令，要使用 str 将其它转换为字符串 | Define the command to run; use str to convert it into a string
    
    TwoPuncture_command         = "./TwoPunctureABE"
    TwoPuncture_command_outfile = "TwoPunctureABE_out.log"

    ## 设置 OpenMP 多线程数目，线程数取输入文件 AMSS_NCKU_Input.py 中设定的 OMP_processes | Set the number of OpenMP threads; the thread count is OMP_processes set in the input file AMSS_NCKU_Input.py
    ## 通过环境变量 OMP_NUM_THREADS 传给 TwoPunctureABE，无需修改 C 程序 | Pass it to TwoPunctureABE through the environment variable OMP_NUM_THREADS, without modifying the C program
    TwoPuncture_env = os.environ.copy()
    TwoPuncture_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    
    TwoPuncture_process = subprocess.Popen(TwoPuncture_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=TwoPuncture_env)
 
    ## 将 TwoPunctureABE 的运行结果写入文件中 | Write the running results of TwoPunctureABE to a file
    
    with open(TwoPuncture_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印 | Read the output in a loop and print it
        for line in TwoPuncture_process.stdout:
            print(line, end='')  # 实时打印输出 | Print the output in real time
            file0.write(line)    # 将行写入文件 | Write the line to the file
            file0.flush()        # 确保每行都被立即写入文件，可选 | Ensure that each line is written to the file immediately (optional)                 
    file0.close()
 
    ## 等待进程结束 | Wait for the process to finish
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

## 这个函数运行 AMSS-NCKU 的后牛顿程序 PNOrbit | This function runs the AMSS-NCKU post-Newtonian program PNOrbit

def run_PNOrbit( input_language ):

    if ( input_language == "Chinese" ):
        print(                                   )
        print( " 正在运行 AMSS-NCKU 程序 PNOrbit " )
        print(                                   ) 
    elif ( input_language == "English" ): 
        print(                                                   ) 
        print( " Running the AMSS-NCKU executable file PNOrbit " ) 
        print(                                                   )
    
    ## 定义要运行的命令，要使用 str 将其它转换为字符串 | Define the command to run; use str to convert it into a string
    
    PNOrbit_command         = "./PNOrbit"
    PNOrbit_command_outfile = "PNOrbit_out.log"

    ## 设置 OpenMP 多线程数目，线程数取输入文件 AMSS_NCKU_Input.py 中设定的 OMP_processes | Set the number of OpenMP threads; the thread count is OMP_processes set in the input file AMSS_NCKU_Input.py
    ## 通过环境变量 OMP_NUM_THREADS 传给 PNOrbit，无需修改 C 程序 | Pass it to PNOrbit through the environment variable OMP_NUM_THREADS, without modifying the C program
    PNOrbit_env = os.environ.copy()
    PNOrbit_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## 使用subprocess.Popen来执行命令，并实时打印输出 | Use subprocess.Popen to execute the command and print the output in real time
    
    PNOrbit_process = subprocess.Popen(PNOrbit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=PNOrbit_env)
 
    ## 将 PNOrbit 的运行结果写入文件中 | Write the running results of PNOrbit to a file
    
    with open(PNOrbit_command_outfile, 'w') as file0:  
        ## 循环读取输出并打印 | Read the output in a loop and print it
        for line in PNOrbit_process.stdout:
            print(line, end='')  # 实时打印输出 | Print the output in real time
            file0.write(line)    # 将行写入文件 | Write the line to the file
            file0.flush()        # 确保每行都被立即写入文件，可选 | Ensure that each line is written to the file immediately (optional)                 
    file0.close()
 
    ## 等待进程结束 | Wait for the process to finish
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
    
