
#################################################
##
## 这个文件包含了数值相对论所需要的输入
## 小曲
## 2024/03/19 --- 2025/09/14
##
#################################################

import numpy   ## 导入 numpy 包 

#################################################

## 设置程序运行目录和计算资源

File_directory   = "BBH_q=1_test_shellpatch"     ## 程序运行目录
Output_directory = "binary_output"               ## 存放二进制数据的子目录
MPI_processes    = 16                            ## 想要调用的进程数目

GPU_Calculation  = "no"                          ## 是否开启 GPU 计算，可选 yes 或 no
CPU_Part         = 0.5
GPU_Part         = 0.5

#################################################


#################################################

## 设置程序计算方法

Symmetry                 = "equatorial-symmetry"   ## 系统对称性，可选 equatorial-symmetry、no-symmetry
Equation_Class           = "BSSN"                  ## 设置方程形式，可选 BSSN、Z4C、BSSN-EScalar、BSSN-EM
                                                   ##       BSSN 和 Z4C   适合于 GR 旋转黑洞的真空计算
                                                   ##       BSSN-EM      涉及 GR 带电黑洞的真空计算
                                                   ##       BSSN-EScalar 涉及到标量张量-F(R) 理论的计算，需要在后面设定额外参数
                                                   ## 注意：GPU 计算仅支持 BSSN
                                                   ## 这里选择 BSSN-EScalar，需要在后面设定 F(R) 理论的参数
Initial_Data_Method      = "Ansorg-TwoPuncture"    ## 设置求解数值相对论初值的方法
                                                   ## 可选 Ansorg-TwoPuncture、
                                                   ##     Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical 
                                                   ## 注意：当前 BSSN-EM 的计算不支持用解析公式 Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical
                                                   ##       当前 BSSN-EScalar 的计算不支持用解析公式 Lousto-Analytical、Cao-Analytical、KerrSchild-Analytical
Time_Evolution_Method    = "runge-kutta-45"        ## 时间演化方法，可选 runge-kutta-45
Finite_Diffenence_Method = "6th-order"             ## 有限差分方法，可选 2nd-order、4th-order、6th-order、8th-order

#################################################


#################################################

## 设置时间演化信息

Start_Evolution_Time     = 0.0                    ## 起始演化时间
Final_Evolution_Time     = 1800.0                 ## 最终演化时间
Check_Time               = 100.0
Dump_Time                = 50.0                   ## 每隔一定时间间隔储存数据
D2_Dump_Time             = 400.0
Analysis_Time            = 0.1
Evolution_Step_Number    = 10000000               ## 时间迭代次数
Courant_Factor           = 0.4                    ## Courant 因子（决定每一步时间演化的时间间隔）
Dissipation              = 0.1                    ## 耗散因子

#################################################


#################################################

## 设置多层格点信息

basic_grid_set    = "Shell-Patch"                    ## 设定网格类型，可选 Patch 和 Shell-Patch
grid_center_set   = "Cell"                           ## 网格中心设置，可选 Cell 和 Vertex

grid_level        = 7                                ## 设置格点的总层数
static_grid_level = 3                                ## 设置静态格点的层数
moving_grid_level = grid_level - static_grid_level   ## 可移动格点的层数

analysis_level    = 0
refinement_level  = 1                                ## 从该层开始进行时间细化

largest_box_xyz_max = [100.0, 100.0, 100.0]          ## 设置最外层格点的坐标最大值
largest_box_xyz_min = - numpy.array(largest_box_xyz_max)  ## 设置最外层格点的坐标最小值

static_grid_number = 96                              ## 设置固定格点每一层每一维数的格点数目（这里对应的 x 轴格点数目，yz 轴格点自动调整）
moving_grid_number = 48                              ## 设置可移动格点每一层每一维数的格点数目
shell_grid_number  = [40, 40, 400]                   ## 设置最外层球状网格（shell patch）的格点数目
                                                     ## 以 phi、theta、r 的顺序给定
devide_factor      = 2.0                             ## 设置相邻两层网格分辨率的比例（不要轻易改变）
static_grid_type   = 'Linear'                        ## 设置固定格点的类型，可选 'Linear'
moving_grid_type   = 'Linear'                        ## 设置固定格点的类型，可选 'Linear'

quarter_sphere_number = 64                           ## 1/4 球面积分的格点数目

#################################################


#################################################

## 设置黑洞 puncture （穿刺法）的信息

puncture_number       = 2                                     ## 设置 puncture 的数目

position_BH           = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的初始位置
parameter_BH          = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的参数
dimensionless_spin_BH = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的无量纲自旋
momentum_BH           = numpy.zeros( (puncture_number, 3) )   ## 初始化每个黑洞的动量

puncture_data_set     = "Manually"                            ## 设置双星轨道坐标的方式，可选 Manually 和 Automatically-BBH

#---------------------------------------------

## 如果设置双星初始轨道坐标的方式选为 Automatically-BBH，只需要给定黑洞参数，偏心率，距离即可

## 这一步与初值求解中的 Ansorg-TwoPuncture 配合使用中需要注意的问题
## 用 Ansorg-TwoPuncture 求解初值，轨道坐标设置可以设置 Manually 和 Automatically-BBH 设置双星轨道坐标
## 但双星轨道坐标如果设置为 Manually 而不是 Automatically-BBH，则要细致设置 Puncture 的位置和动量取值，否则可能会使 TwoPuncture 程序无法正确读入输入而报错）

Distance = 11.0
e0       = 0.0

## 设置每个黑洞的参数 (M Q* a*)  
## 质量  无量纲电荷  无量纲自旋
parameter_BH[0] = [ 0.487208758,  0.0,  0.0 ]   
parameter_BH[1] = [ 0.487208758,  0.0, -0.0 ]  
## 注意，如果求解数值相对论初值的方法选为 Ansorg-TwoPuncture，第一个黑洞必须为质量较大的那个，且黑洞总质量会自动 rescale 为 M=1 （其它情况下必须手动 rescale）

## 设置每个黑洞的无量纲自旋
## 无对称性时 ，需要手动给 3 个方向的自旋角动量
dimensionless_spin_BH[0] = [ 0.0,  0.0,  0.0 ]   
dimensionless_spin_BH[1] = [ 0.0,  0.0, -0.0 ]  

## 注意，如果设置双星初始轨道坐标的方式选为 Automatically-BBH，则程序自动调整将较大质量黑洞放在 y 轴正向，将较小质量黑洞放在 y 轴负向
##       如果设置双星初始轨道坐标的方式选为 Manually，则需要手动调整到 y 轴方向 
## use Brugmann's convention
##  -----0-----> y
##   -      +     

#---------------------------------------------

## 如果设置 puncture 初始轨道坐标的方式选为 Manually，还需要手动给定所有黑洞参数

## 设置每个黑洞的初始位置
position_BH[0]  = [  0.0,  +5.5,  0.0 ]  
position_BH[1]  = [  0.0,  -5.5,  0.0 ]  

## 设置每个黑洞的动量信息  
momentum_BH[0]  = [ -0.090109887, -0.000703975,   0.0 ]
momentum_BH[1]  = [ +0.090109887, +0.000703975,   0.0 ] 


#################################################


#################################################

## 设置引力波和探测器的相关信息

GW_L_max        = 4                      ## 引力波最大的 L
GW_M_max        = 4                      ## 引力波最大的 M
Detector_Number = 11                     ## 探测器的数目
Detector_Rmin   = 50.0                   ## 最近探测器的距离
Detector_Rmax   = 150.0                  ## 最远探测器的距离

#################################################


#################################################

## 设置表观视界的参数

AHF_Find       = "no"                    ## 是否开启表观视界计算，可选 yes 或 no

AHF_Find_Every = 24
AHF_Dump_Time  = 20.0

#################################################


#################################################

## 标量-张量-f(R) 理论的一些参数
## 仅对 BSSN-EScalar 的计算有影响

FR_a2     = 3.0        ## f(R) = R + a2 * R^2    
FR_l2     = 10000.0
FR_phi0   = 0.00005
FR_r0     = 120.0
FR_sigma0 = 8.0
FR_Choice = 2          ## Choice 可选为 1 2 3 4 5
                       ## 1: phi(r) = phi0 * Exp(-(r-r0)**2/sigma0)   
                       ##    V(r)   = 0
                       ## 2: phi(r) =  phi0 * a2^2/(1+a2^2)  
                       ##    V(r)   = Exp(-8*Sqrt(PI/3)*phi(r)) * (1-Exp(4*Sqrt(PI/3)*phi(r)))**2 / (32*PI*a2)
                       ##    该 V(r) 由  f(R) = R + a2*R^2 诱导
                       ## 3: Schrodinger-Newton 系统给定的 phi(r) 
                       ##    V(r)   = Exp(-8*Sqrt(PI/3)*phi(r)) * (1-Exp(4*Sqrt(PI/3)*phi(r)))**2 / (32*PI*a2)
                       ##    该 V(r) 由  f(R) = R + a2*R^2 诱导
                       ## 4: phi(r) = phi0 * 0.5 * ( tanh((r+r0)/sigma0) - tanh((r-r0)/sigma0) )  
                       ##    V(r)   = 0
                       ##    f(R)   = R + a2*R^2  其中 a2 设定为 a2 = +oo
                       ## 5: phi(r) = phi0 * Exp(-(r-r0)**2/sigma)   
                       ##    V(r)   = 0

#################################################


#################################################

## 其它选项
## 还在测试中
## 但不建议用户轻易改动这些选项

boundary_choice = "BAM-choice"     ## 索莫菲边界条件设定，可选 "BAM-choice" 和 "Shibata-choice"
                                   ## 目前的版本定建议选为 "BAM-choice"

gauge_choice  = 2                  ## 规范条件选取
                                   ## 0: B^i gauge
                                   ## 1: David's puncture gauge
                                   ## 2: MB B^i gauge               ## 对Z4C和GPU计算好像有bug
                                   ## 3: RIT B^i gauge
                                   ## 4: MB beta gauge 
                                   ## 5: RIT beta gauge 
                                   ## 6: MGB1 B^i gauge
                                   ## 7: MGB2 B^i gauge
                                   ## 目前的版本建议选为 0 或 1
                                   
tetrad_type  = 2                   ## tetradtype 选取
                                   ## 以下   v:r; u: phi; w: theta
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
                                   ## 目前的版本建议选为 2
                                   
#################################################
                                   
