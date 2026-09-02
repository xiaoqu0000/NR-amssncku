
##################################################################
##
## 这个文件设定了 AMSS-NCKU 编译和运行的相关命令
## 小曲
## 2025/01/24
## 2026/06/01 参考浙江大学CPU优化版设定
##
##################################################################


import os
import shutil
import subprocess

import AMSS_NCKU_Input as input_data


def _resolve_mpi_launcher():
    """Return (absolute_mpirun_path, flavor) where flavor is 'intel' or 'generic'.

    Resolution order:
      1. $MPI_LAUNCHER env var if set to a non-'auto' value.
      2. `mpirun` on PATH.
    Raises RuntimeError if no launcher can be found.
    """
    candidate = os.environ.get("MPI_LAUNCHER", "auto")
    if candidate and candidate != "auto":
        path = candidate
        if not os.path.isabs(path):
            resolved = shutil.which(path)
            if resolved:
                path = resolved
    else:
        path = shutil.which("mpirun")
    if not path or not os.path.exists(path):
        raise RuntimeError(
            "Could not locate an mpirun binary. Install OpenMPI/MPICH "
            "(e.g. `sudo apt-get install openmpi-bin`) or set MPI_LAUNCHER "
            "in config.env to the absolute path of your launcher."
        )
    flavor = "generic"
    try:
        version = subprocess.check_output(
            [path, "--version"], stderr=subprocess.STDOUT, text=True, timeout=10
        )
    except (subprocess.CalledProcessError, subprocess.TimeoutExpired, OSError):
        version = ""
    if "Intel" in version and "MPI" in version:
        flavor = "intel"
    return path, flavor


##################################################################



##################################################################

## Generate processor binding list for Intel MPI

def generate_pin_processor_list(mpi_processes, omp_processes, cores_per_numa):
    """
    生成 Intel MPI 的 I_MPI_PIN_PROCESSOR_LIST 绑核字符串
    
    参数:
        mpi_processes: MPI 进程数
        omp_processes: 每个 MPI 进程的 OMP 线程数
        cores_per_numa: 单个 NUMA 节点的核数
    
    返回:
        绑核字符串：
        - 单 MPI 进程: "0,1,2,3,4,5,6,7"
        - 多 MPI 进程: "'(0,1,2,3,4,5,6,7),(8,9,10,11,12,13,14,15)'"
        
    示例:
        generate_pin_processor_list(1, 8, 8)  -> "0,1,2,3,4,5,6,7"
        generate_pin_processor_list(2, 8, 8)  -> "'(0,1,2,3,4,5,6,7),(8,9,10,11,12,13,14,15)'"
        generate_pin_processor_list(4, 4, 8)  -> "'(0,1,2,3),(8,9,10,11),(16,17,18,19),(24,25,26,27)'"
    """
    if mpi_processes == 1:
        # 单 MPI 进程：不加括号，直接返回核心列表
        start_core = 0
        core_list = [str(start_core + i) for i in range(omp_processes)]
        return ','.join(core_list)
    else:
        # 多 MPI 进程：每组加括号，外面再加单引号
        processor_groups = []
        
        for mpi_rank in range(mpi_processes):
            # 每个 MPI 进程放在不同的 NUMA 节点上
            numa_id = mpi_rank
            start_core = numa_id * cores_per_numa
            
            # 为该 MPI 进程分配的核心列表
            core_list = [str(start_core + i) for i in range(omp_processes)]
            processor_groups.append(f"({','.join(core_list)})")
        
        return ','.join(processor_groups)

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

    ## Build command
    ## Timing flag: pass EXTRA_CXXFLAGS to make so it gets appended to CXXAPPFLAGS
    timing_flag = "-DBSSN_TIMING_ENABLED" if getattr(input_data, "Enable_Timing", "no") == "yes" else ""
    if (input_data.GPU_Calculation == "no"):
        makefile_command  = f"make -j8 ABE EXTRA_CXXFLAGS='{timing_flag}' EXTRA_F90FLAGS='{timing_flag}'"
    elif (input_data.GPU_Calculation == "yes"):
        makefile_command  = f"make -j8 ABEGPU EXTRA_CXXFLAGS='{timing_flag}' EXTRA_F90FLAGS='{timing_flag}'"
    else:
        if ( input_language == "Chinese" ):
            print( " CPU/GPU 计算设置出错 " )
            print(                         )
        elif ( input_language == "English" ):
            print( " CPU/GPU numerical calculation setting is wrong " )
            print(                                                    )
 
    ## Execute the command with subprocess.Popen and stream output
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    ## Read and print output lines as they arrive
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出

    ## Wait for the process to finish
    makefile_return_code = makefile_process.wait()
    if makefile_return_code != 0:
        raise subprocess.CalledProcessError(makefile_return_code, makefile_command)
        
    if ( input_language == "Chinese" ):
        print(                               )
        print( " AMSS-NCKU 程序 ABE 编译完成" ) 
        print(                               )
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
    
    ## Build command
    makefile_command = "make" + " TwoPunctureABE"

    ## Execute the command with subprocess.Popen and stream output
    makefile_process = subprocess.Popen(makefile_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True) 
    
    ## Read and print output lines as they arrive
    for line in makefile_process.stdout:
        print(line, end='')  # 实时打印输出
        
    ## Wait for the process to finish
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
        print(                                   )
        print( " 正在编译 AMSS-NCKU 程序 PNOrbit " ) 
        print(                                   )
    elif ( input_language == "English" ):
        print(                                                     )
        print( " Compiling the AMSS-NCKU executable file PNOrbit " )
        print(                                                     )
    
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
        print(                                  )
        print( " AMSS-NCKU 程序 PNOrbit 编译完成" ) 
        print(                                  )
    elif ( input_language == "English" ):   
        print(                                                                      ) 
        print( " Compilation of the AMSS-NCKU executable file PNOrbit is finished " )
        print(                                                                      )
    
    return
    
##################################################################



##################################################################

## Run the AMSS-NCKU main program ABE

def run_ABE( input_language ):

    if ( input_language == "Chinese" ):
        print(                                        )
        print( " 正在运行 AMSS-NCKU 主程序 ABE/ABEGPU " ) 
        print(                                        )
    elif ( input_language == "English" ): 
        print(                                                      )    
        print( " Running the AMSS-NCKU executable file ABE/ABEGPU " ) 
        print(                                                      )

    ## 定义要运行的命令，要将其转换为字符串
    
    if (input_data.GPU_Calculation == "no"):
        # Resolve mpirun at runtime instead of hardcoding the Intel MPI path,
        # so the same script works on plain OpenMPI installs and on Intel MPI.
        mpirun_path, mpi_flavor = _resolve_mpi_launcher()

        pin_list = generate_pin_processor_list(
            # input_data.MPI_processes,
            1,
            input_data.OMP_processes * input_data.MPI_processes,
            input_data.cores_per_numa * input_data.MPI_processes
        )

        print(f" MPI launcher: {mpirun_path} (flavor: {mpi_flavor})")
        print(f" Processor binding list: {pin_list}")
        print()

        # OMP env vars are portable across MPI implementations.
        # OMP_DISPLAY_ENV and OMP_DISPLAY_AFFINITY are explicitly set to FALSE
        # (overriding any exported value) so that the OpenMP runtime does NOT
        # print the environment / thread affinity info at program startup.
        omp_env = (
            f"-x OMP_NUM_THREADS={input_data.OMP_processes} "
            f"-x OMP_PROC_BIND=close "
            f"-x OMP_PLACES=cores "
            f"-x OMP_WAIT_POLICY=ACTIVE "
            f"-x GOMP_SPINCOUNT=INFINITE "
            f"-x OMP_SCHEDULE=static "
            f"-x OMP_DISPLAY_ENV=FALSE "
            f"-x OMP_DISPLAY_AFFINITY=FALSE "
        )

        if mpi_flavor == "intel":
            # Intel MPI uses `-genv` and understands I_MPI_PIN_* knobs.
            intel_env = omp_env.replace("-x ", "-genv ")
            mpi_command = (
                f"{mpirun_path} "
                f"{intel_env}"
                f"-genv I_MPI_PIN_DOMAIN=omp "
                f"-genv I_MPI_PIN_PROCESSOR_LIST={pin_list} "
                f"-np {input_data.MPI_processes} ./ABE"
            )
        else:
            # OpenMPI / MPICH: use `-x` for env propagation and standard binding flags.
            mpi_command = (
                f"{mpirun_path} "
                f"{omp_env}"
                f"--bind-to core --map-by slot:pe={input_data.OMP_processes} "
                f"-np {input_data.MPI_processes} ./ABE"
            )
        mpi_command_outfile = "ABE_out.log"
    elif (input_data.GPU_Calculation == "yes"):
        mpi_command         = "mpirun -np " + str(input_data.MPI_processes) + " ./ABEGPU"
        mpi_command_outfile = "ABEGPU_out.log"
 
    print(mpi_command)

    ## Execute the MPI command and stream output
    mpi_process = subprocess.Popen(mpi_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)

    ## Write ABE run output to file while printing to stdout
    with open(mpi_command_outfile, 'w') as file0:  
        ## Read and print output lines; also write each line to file
        for line in mpi_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # write the line to file
            file0.flush()        # flush to ensure each line is written immediately (optional)            
    file0.close()

    ## Wait for the process to finish
    mpi_return_code = mpi_process.wait()
    
    if ( input_language == "Chinese" ):
        print(                                         )
        print( " AMSS-NCKU 主程序 ABE/ABEGPU 运行结束 " ) 
        print(                                         )
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
    
    ## Define the command to run
    TwoPuncture_command         = "./TwoPunctureABE"
    TwoPuncture_command_outfile = "TwoPunctureABE_out.log"

    ## Set the number of OpenMP threads from OMP_processes in AMSS_NCKU_Input.py
    ## via the environment variable OMP_NUM_THREADS (no C code modification needed)
    TwoPuncture_env = os.environ.copy()
    TwoPuncture_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## Execute the command with subprocess.Popen and stream output
    TwoPuncture_process = subprocess.Popen(TwoPuncture_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=TwoPuncture_env)

    ## Write TwoPunctureABE run output to file while printing to stdout
    with open(TwoPuncture_command_outfile, 'w') as file0:  
        ## Read and print output lines; also write each line to file
        for line in TwoPuncture_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # write the line to file
            file0.flush()        # flush to ensure each line is written immediately (optional)                 
    file0.close()

    ## Wait for the process to finish
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
    
    ## Define the command to run
    PNOrbit_command         = "./PNOrbit"
    PNOrbit_command_outfile = "PNOrbit_out.log"

    ## Set the number of OpenMP threads from OMP_processes in AMSS_NCKU_Input.py
    ## via the environment variable OMP_NUM_THREADS (no C code modification needed)
    PNOrbit_env = os.environ.copy()
    PNOrbit_env['OMP_NUM_THREADS'] = str(input_data.OMP_processes)

    ## Execute the command with subprocess.Popen and stream output
    PNOrbit_process = subprocess.Popen(PNOrbit_command, shell=True, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=PNOrbit_env)

    ## Write PNOrbit run output to file while printing to stdout
    with open(PNOrbit_command_outfile, 'w') as file0:  
        ## Read and print output lines; also write each line to file
        for line in PNOrbit_process.stdout:
            print(line, end='')  # 实时打印输出
            file0.write(line)    # write the line to file
            file0.flush()        # flush to ensure each line is written immediately (optional)                 
    file0.close()

    ## Wait for the process to finish
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
    
