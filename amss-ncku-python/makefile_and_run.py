
##################################################################
##
## 这个文件设定了 AMSS-NCKU 编译和运行的相关命令
## 小曲
## 2025/01/24
##
##################################################################


import AMSS_NCKU_Input as input_data
import subprocess


##################################################################



##################################################################

## 这个函数编译 AMSS-NCKU 主程序 ABE

def makefile_ABE():

    print(                                                        )
    print( " 正在编译 AMSS-NCKU 程序 ABE/ABEGPU "                   )   
    print( " Compiling the AMSS-NCKU executable file ABE/ABEGPU " ) 
    print(                                                        )

    ## 编译命令
    if (input_data.GPU_Calculation == "no"):
        makefile_command  = "make -j4" + " ABE"
    elif (input_data.GPU_Calculation == "yes"):
        makefile_command  = "make -j4" + " ABEGPU"
    else:
        print( " CPU/GPU 计算设置出错 "                             )
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
        
    print(                                                                  )
    print( " AMSS-NCKU 程序 ABE 编译完成"                                     ) 
    print( " Compilation of the AMSS-NCKU executable file ABE is finished " ) 
    print(                                                                  )
    
    return
        
##################################################################



##################################################################

## 这个函数编译 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE

def makefile_TwoPunctureABE():

    print(                                                            )
    print( " 正在编译 AMSS-NCKU 程序 TwoPunctureABE "                   ) 
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
        
    print(                                                                             )
    print( " AMSS-NCKU 程序 TwoPunctureABE 编译完成"                                     )     
    print( " Compilation of the AMSS-NCKU executable file TwoPunctureABE is finished " )
    print(                                                                             )
    
    return
    
##################################################################



##################################################################

## 这个函数运行 AMSS-NCKU 主程序 ABE

def run_ABE():

    print(                                                      )
    print( " 正在运行 AMSS-NCKU 主程序 ABE/ABEGPU "               )      
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
    
    print(                                           )
    print( " AMSS-NCKU 主程序 ABE/ABEGPU 运行结束 "    ) 
    print( " The ABE/ABEGPU simulation is finished " ) 
    print(                                           )
    
    return

##################################################################



##################################################################

## 这个函数运行 AMSS-NCKU 的 TwoPuncture 程序 TwoPunctureABE

def run_TwoPunctureABE():

    print(                                                          )
    print( " 正在运行 AMSS-NCKU 程序 TwoPunctureABE "                 )   
    print( " Running the AMSS-NCKU executable file TwoPunctureABE " ) 
    print(                                                          )
    
    ## 定义要运行的命令，要使用 str 将其它转换为字符串
    
    TwoPuncture_command         = "./TwoPunctureABE"
    TwoPuncture_command_outfile = "TwoPunctureABE_out.log"

    ## 使用subprocess.Popen来执行命令，并实时打印输出
    
    TwoPuncture_process = subprocess.Popen(TwoPuncture_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
 
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
    
    print(                                               )
    print( " AMSS-NCKU 程序 TwoPunctureABE 运行结束 "      ) 
    print( " The TwoPunctureABE simulation is finished " ) 
    print(                                               )
    
    return

##################################################################
    
