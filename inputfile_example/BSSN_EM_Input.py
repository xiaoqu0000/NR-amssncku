
#################################################
##
## 这个文件包含了数值相对论所需要的输入 | This file contains the inputs required for numerical relativity
## 小曲 | Xiao Qu
## 2024/03/19 --- 2026/09/25
##
#################################################

import numpy   

#################################################

## 设置程序运行目录和计算资源 | Set the program running directory and computing resources

File_directory   = "GW150914_testBSSNEM"         ## 程序运行目录 | Program running directory
Output_directory = "output_file"                 ## 存放二进制数据的子目录 | Subdirectory for storing binary data

MPI_processes    = 16                            ## 想要调用的进程数目        number of mpi processes used in the simulation
OMP_processes    = 4                             ## 想要调用的线程数目        number of OpenMP threads per MPI process

cores_per_numa   = OMP_processes                 ## number of CPU cores per NUMA node (for processor binding)

GPU_Calculation  = "no"                          ## 是否开启 GPU 计算，可选 yes 或 no   Use GPU or not 
GPU_Code_Choice  = "flash"                       ## GPU 内层代码选择，可选 flash 和 pro | GPU inner code choice, options: flash and pro 

#################################################


#################################################

## 设置程序计算方法 | Set the computation methods of the program

Symmetry                 = "equatorial-symmetry"   ## 系统对称性，可选 equatorial-symmetry、no-symmetry | System symmetry, options: equatorial-symmetry, no-symmetry
Equation_Class           = "BSSN-EM"               ## 设置方程形式，可选 BSSN、Z4C、BSSN-EScalar、BSSN-EM | Set the form of the equations, options: BSSN, Z4C, BSSN-EScalar, BSSN-EM
                                                   ##       BSSN 和 Z4C   适合于 GR 旋转黑洞的真空计算 | suitable for vacuum computations of rotating black holes in GR
                                                   ##       BSSN-EM      涉及 GR 带电黑洞的真空计算 | involves vacuum computations of charged black holes in GR
                                                   ##       BSSN-EScalar 涉及到标量张量-F(R) 理论的计算，需要在后面设定额外参数 | involves computations in the scalar-tensor f(R) theory, additional parameters need to be set below
                                                   ## 这里没有选择 BSSN-EScalar， F(R) 理论的参数不是必须要设定的 | BSSN-EScalar is not chosen here, so the parameters of the F(R) theory are not mandatory
Initial_Data_Method      = "Ansorg-TwoPuncture"    ## 设置求解数值相对论初值的方法 | Set the method for solving the initial data of numerical relativity
                                                   ## 可选 Ansorg-TwoPuncture、Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical | Options: Ansorg-TwoPuncture, Lousto-Analytical, Cao-Analytical, KerrSchild-Analytical
                                                   ## 注意：当前 BSSN-EM 的计算不支持用解析公式 Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical | Note: the current BSSN-EM computation does not support the analytical formulas Lousto-Analytical, Cao-Analytical, KerrSchild-Analytical
                                                   ##       当前 BSSN-EScalar 的计算不支持用解析公式 Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical | the current BSSN-EScalar computation does not support the analytical formulas Lousto-Analytical, Cao-Analytical, KerrSchild-Analytical
Time_Evolution_Method    = "runge-kutta-45"        ## 时间演化方法，可选 runge-kutta-45 | Time evolution method, options: runge-kutta-45
Finite_Diffenence_Method = "6th-order"             ## 有限差分方法，可选 2nd-order、4th-order、6th-order、8th-order | Finite-difference method, options: 2nd-order, 4th-order, 6th-order, 8th-order

#################################################


#################################################

## 设置时间演化信息 | Set the time evolution information

Start_Evolution_Time     = 0.0                    ## 起始演化时间 | Start time of the evolution
Final_Evolution_Time     = 1500.0                 ## 最终演化时间 | Final time of the evolution
Check_Time               = 100.0
Dump_Time                = 50.0                   ## 每隔一定时间间隔储存数据 | Dump data at regular time intervals
D2_Dump_Time             = 300.0
Analysis_Time            = 0.1
Evolution_Step_Number    = 10000000               ## 时间迭代次数 | Number of time iteration steps
Courant_Factor           = 0.4                    ## Courant 因子（决定每一步时间演化的时间间隔） | Courant factor (determines the time interval of each evolution step)
Dissipation              = 0.2                    ## 耗散因子 | Dissipation factor

#################################################


#################################################

## 设置多层格点信息 | Set the multilayer grid information

basic_grid_set    = "Patch"                          ## 设定网格类型，可选 Patch 和 Shell-Patch | Set the grid type, options: Patch and Shell-Patch
grid_center_set   = "Cell"                           ## 网格中心设置，可选 Cell 和 Vertex | Grid center setting, options: Cell and Vertex

grid_level        = 10                               ## 设置格点的总层数 | Total number of grid levels
static_grid_level = 6                                ## 设置静态格点的层数 | Number of static grid levels
moving_grid_level = grid_level - static_grid_level   ## 可移动格点的层数 | Number of moving grid levels

analysis_level    = 0
refinement_level  = 4                                ## 从该层开始进行时间细化 | Time refinement starts from this level

largest_box_xyz_max = [500.0, 500.0, 500.0]          ## 设置最外层格点的坐标最大值 | Maximum coordinates of the outermost grid level
largest_box_xyz_min = - numpy.array(largest_box_xyz_max)  ## 设置最外层格点的坐标最小值 | Minimum coordinates of the outermost grid level

static_grid_number = 128                             ## 设置固定格点每一层每一维数的格点数目（这里对应的 x 轴格点数目，yz 轴格点自动调整） | Number of grid points per dimension per static level (given here for the x axis; grid points in y and z are adjusted automatically)
moving_grid_number = 64                              ## 设置可移动格点每一层每一维数的格点数目 | Number of grid points per dimension per moving level
shell_grid_number  = [32, 32, 100]                   ## 设置最外层球状网格（shell patch）的格点数目 | Number of grid points of the outermost spherical grid (shell patch)
                                                     ## 以 phi、theta、r 的顺序给定 | Given in the order phi, theta, r
devide_factor      = 2.0                             ## 设置相邻两层网格分辨率的比例（不要轻易改变） | Resolution ratio between adjacent grid levels (do not change it lightly)
static_grid_type   = 'Linear'                        ## 设置固定格点的类型，可选 'Linear' | Type of the static grid, options: 'Linear'
moving_grid_type   = 'Linear'                        ## 设置固定格点的类型，可选 'Linear' | Type of the static grid, options: 'Linear'

quarter_sphere_number = 128                          ## 1/4 球面积分的格点数目 | Number of grid points for the 1/4-sphere surface integral

#################################################


#################################################

## 设置黑洞 puncture （穿刺法）的信息 | Set the black hole puncture (puncture method) information

puncture_number       = 2                                     ## 设置 puncture 的数目 | Number of punctures

position_BH           = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的初始位置 | Initialize the initial position of each black hole
parameter_BH          = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的参数 | Initialize the parameters of each black hole
dimensionless_spin_BH = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的无量纲自旋 | Initialize the dimensionless spin of each black hole
momentum_BH           = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的动量 | Initialize the momentum of each black hole

puncture_data_set     = "Manually"                            ## 设置双星轨道坐标的方式，可选 Manually 和 Automatically-BBH | Method for setting the binary orbital coordinates, options: Manually and Automatically-BBH

##---------------------------------------------
## 如果设置双星初始轨道坐标的方式选为 Automatically-BBH，只需要给定黑洞参数，偏心率，距离即可 | If the binary initial orbital coordinates are set via Automatically-BBH, only the black hole parameters, eccentricity, and separation are needed
##---------------------------------------------

## 这一步与初值求解中的 Ansorg-TwoPuncture 配合使用中需要注意的问题 | Points to note when using this step together with Ansorg-TwoPuncture in the initial data solve
## 用 Ansorg-TwoPuncture 求解初值，轨道坐标设置可以设置 Manually 和 Automatically-BBH 设置双星轨道坐标 | When solving initial data with Ansorg-TwoPuncture, the binary orbital coordinates can be set via either Manually or Automatically-BBH
## 但双星轨道坐标如果设置为 Manually 而不是 Automatically-BBH，则要细致设置 Puncture 的位置和动量取值，否则可能会使 TwoPuncture 程序无法正确读入输入而报错） | If the binary orbital coordinates are set to Manually rather than Automatically-BBH, the puncture positions and momenta must be set carefully, otherwise the TwoPuncture code may fail to read the input correctly and abort with an error)

## 双黑洞系统的初始间距
## initial orbital distance and ellipticity for BBHs system
## ( needed for "Automatically-BBH" case , not affect the "Manually" case )
Allow_PN_Evaluation = "no"  ## 是否允许从更远的位置进行后牛顿演化，来设定初值 | Whether to allow a post-Newtonian evolution from a farther position to set the initial data
Distance_initial = 60.0
e0               = 0.0
Distance_final   = 10.0

## 设置每个黑洞的参数 (M Q* a*) | Set the parameters of each black hole (M Q* a*)  
## 质量  无量纲电荷  无量纲自旋 | Mass, dimensionless charge, dimensionless spin
parameter_BH[0] = [ 36.0/(36.0+29.0),  0.5,  0.31 ]   
parameter_BH[1] = [ 29.0/(36.0+29.0),  0.5, -0.46 ]  
## 注意，如果求解数值相对论初值的方法选为 Ansorg-TwoPuncture，第一个黑洞必须为质量较大的那个，且黑洞总质量会自动 rescale 为 M=1 （其它情况下必须手动 rescale） | Note: if the initial data method is Ansorg-TwoPuncture, the first black hole must be the more massive one, and the total mass is automatically rescaled to M=1 (in other cases the rescaling must be done manually)

## 设置每个黑洞的无量纲自旋 | Set the dimensionless spin of each black hole
## 无对称性时 ，需要手动给 3 个方向的自旋角动量 | For no symmetry, the spin angular momenta in all 3 directions must be given manually
dimensionless_spin_BH[0] = [ 0.0,  0.0,  0.31 ]   
dimensionless_spin_BH[1] = [ 0.0,  0.0, -0.46 ]  

## 注意，如果设置双星初始轨道坐标的方式选为 Automatically-BBH，则程序自动调整将较大质量黑洞放在 y 轴正向，将较小质量黑洞放在 y 轴负向 | Note: if the binary initial orbital coordinates are set via Automatically-BBH, the code automatically places the more massive black hole on the positive y axis and the less massive one on the negative y axis
##       如果设置双星初始轨道坐标的方式选为 Manually，则需要手动调整到 y 轴方向 | if the binary initial orbital coordinates are set via Manually, the coordinates must be adjusted manually to the y axis direction 
## use Brugmann's convention
##  -----0-----> y
##   -      +     

##---------------------------------------------
## 如果设置 puncture 初始轨道坐标的方式选为 Manually，还需要手动给定所有黑洞参数 | If the puncture initial orbital coordinates are set via Manually, all black hole parameters must still be given manually
##---------------------------------------------

## 设置每个黑洞的初始位置 | Set the initial position of each black hole
position_BH[0]  = [  0.0,  +4.4615385,  0.0 ]  
position_BH[1]  = [  0.0,  -5.5384615,  0.0 ]  

## 设置每个黑洞的动量信息 | Set the momentum information of each black hole  
momentum_BH[0]  = [ -0.0953015, -0.00084515,   0.0 ]
momentum_BH[1]  = [ +0.0953015, +0.00084515,   0.0 ]


#################################################


#################################################

## 设置引力波和探测器的相关信息 | Set the gravitational wave and detector information

GW_L_max        = 4                      ## 引力波最大的 L | Maximum L of the gravitational wave
GW_M_max        = 4                      ## 引力波最大的 M | Maximum M of the gravitational wave
Detector_Number = 11                     ## 探测器的数目 | Number of detectors
Detector_Rmin   = 50.0                   ## 最近探测器的距离 | Distance of the nearest detector
Detector_Rmax   = 150.0                  ## 最远探测器的距离 | Distance of the farthest detector

#################################################


#################################################

## 设置表观视界的参数 | Set the apparent horizon parameters

AHF_Find       = "no"                    ## 是否开启表观视界计算，可选 yes 或 no | Whether to enable apparent horizon finding, options: yes or no

AHF_Find_Every = 24
AHF_Dump_Time  = 20.0

#################################################


#################################################

## 设置二进制画图相关信息 | Set the binary data plotting information

plot_binary_data_set         = "xy-xz-yz-plot"             ## 选择 "xy-plot", "xz-plot", "yz-plot", "xy-xz-yz-plot" | choose one of "xy-plot", "xz-plot", "yz-plot", "xy-xz-yz-plot"
## plot_binary_data_type        = "Linear"                 ## 对哈密顿约束，TrK这些数据画图选择 | Plotting options for data such as the Hamiltonian constraint and TrK
                                                           ## 选择 "Linear", "Logarithmic" | choose one of "Linear", "Logarithmic"
                                                           ## 后续开发 | To be developed later
plot_binary_data_level       = "All-Level"                 ## 选择 "All-level", "Single-level" | choose one of "All-level", "Single-level"
plot_binary_data_levelnumber = static_grid_level - 1       ## 选择对哪个层来进行画图 | Choose which level to plot
plot_binary_data_processes   = OMP_processes               ## 二进制数据并行画图使用的进程数目 | number of processes used in parallel binary data plotting
                                                           ## 如果设为 0，则自动使用所有 CPU 核数 | if set to 0, all CPU cores are used automatically


#################################################


#################################################

## 其它选项 | Other options
## 还在测试中 | Still under testing
## 但不建议用户轻易改动这些选项 | But it is not recommended to change these options lightly

boundary_choice = "BAM-choice"     ## 索莫菲边界条件设定，可选 "BAM-choice" 和 "Shibata-choice" | Sommerfeld boundary condition setting, options: "BAM-choice" and "Shibata-choice"
                                   ## 目前的版本定建议选为 "BAM-choice" | The current version recommends "BAM-choice"

gauge_choice  = 0                  ## 规范条件选取 | Gauge condition choice
                                   ## 0: B^i gauge
                                   ## 1: David's puncture gauge
                                   ## 2: MB B^i gauge               ## 对Z4C和GPU计算好像有bug | Seems to have a bug for Z4C and GPU computations
                                   ## 3: RIT B^i gauge
                                   ## 4: MB beta gauge 
                                   ## 5: RIT beta gauge 
                                   ## 6: MGB1 B^i gauge
                                   ## 7: MGB2 B^i gauge
                                   ## 目前的版本建议选为 0 或 1 | The current version recommends 0 or 1
                                   
tetrad_type  = 2                   ## tetradtype 选取 | Tetrad type choice
                                   ## 以下   v:r; u: phi; w: theta | below, v: r; u: phi; w: theta
                                   ##      v^a = (x,y,z)
                                   ## 0: orthonormal order: v,u,w
                                   ##    v^a = (x,y,z)   
                                   ##    m = (phi - i theta)/sqrt(2) 
                                   ##    following Frans, Eq.(8) of  PRD 75, 124018(2007)
                                   ## 1: orthonormal order: w,u,v
                                   ##    m = (theta + i phi)/sqrt(2) 
                                   ##    following Sperhake, Eq.(3.2) of  PRD 85, 124062(2012)    
                                   ## 2: orthonormal order: v,u,w
                                   ##    v_a = (x,y,z)
                                   ##    m = (phi - i theta)/sqrt(2) 
                                   ##    following Frans, Eq.(8) of  PRD 75, 124018(2007)
                                   ## 目前的版本建议选为 2 | The current version recommends 2
                                   
#################################################

                                   
