
#################################################
##
## 这个文件对 AMSS-NCKU 数值相对论的结果中提取引力波振幅，并进行画图
## 小曲
## 2024/10/01 --- 2025/11/20
##
#################################################

import numpy                               ## 导入 numpy 包进行数组的操作
import scipy                               ## 导入进行插值
import math
import matplotlib.pyplot    as     plt     ## 导入 matplotlib 包进行画图
import os                                  ## 导入 os 包进行系统操作

import AMSS_NCKU_Input as input_data

# plt.rcParams['text.usetex'] = True  ## 在绘图中允许使用 latex 字体


####################################################################################



####################################################################################

## 将 [t0, t1] 范围的离散时间数据进行傅里叶变换
## 根据 deepseek 例子修改

## 计算非周期离散数据的傅里叶变换
    
## 参数：
## time : 时间序列 [0, t0]
## f :    函数值序列
## apply_window : 是否加窗（推荐True）
## zero_pad_factor : 零填充倍数（推荐2-8）
    
## 返回：
## frequency : 频率轴（Hz）
## frequency_spectrum : 复数频谱

def compute_frequency_spectrum(time, signal_f, apply_window=True, zero_pad_factor=4):

    ## 计算采样率
    N  = len(time)
    dt = time[1] - time[0]        ## 采样间隔
    fs = 1/dt                     ## 采样频率
    omega_s = fs * 2.0 * math.pi  ## 采样角频率
    
    ## 数据预处理
    ## f_detrended = signal_f - numpy.mean(signal_f)  ## 去直流分量
    f_detrended = signal_f
    
    ## 加窗处理（抑制频谱泄漏）
    if apply_window:
        # window = scipy.signal.windows.hann(N)  # 汉宁窗
        window = scipy.signal.windows.tukey(N, alpha=0.1)  # 或Tukey窗
        f_windowed = f_detrended * window
    else:
        f_windowed = f_detrended
    
    # 零填充（提高频率分辨率）
    M = zero_pad_factor * N
    f_padded = numpy.zeros(M)
    f_padded[:N] = f_windowed
    
    # 执行傅里叶变换 FFT，得到复数频谱

    # 这里用 numpy.fft.fft 进行快速傅里叶变换
    # 基本语法为 numpy.fft.fft(a, n=None, axis=-1, norm=None)
    # a：输入的一维数组，代表要进行傅里叶变换的信号。
    # n：可选参数，指定傅里叶变换的长度。如果 n 大于 a 的长度，会在 a 的末尾补零；如果 n 小于 a 的长度，则会截断 a。默认值为 a 的长度。
    # axis：可选参数，指定在哪个轴上进行傅里叶变换，默认是最后一个轴。
    # norm：可选参数，指定归一化方式，取值可以是 None、'ortho' 等
    # numpy.fft.fft 函数返回的结果是一个复数数组，其频率轴的范围是从 -fs/2 到 fs/2（fs 是采样率）
    
    frequency_spectrum = numpy.fft.fft(f_padded, norm='ortho')    
    # frequency_spectrum = numpy.fft.fft(f_detrended, norm='ortho')  
    # frequency_spectrum = numpy.fft.fft(f_windowed, norm='ortho')     
    
    # 频率轴生成

    # 这里用 numpy.fft.fftfreq 生成离散傅里叶频谱的频率轴
    # 其函数签名为 numpy.fft.fftfreq(n, d=1.0)，参数含义如下：
    # n：这是一个必需的整数参数，表示输入信号的长度，也就是进行傅里叶变换时输入序列的样本点数
    # d：这是一个可选的浮点数参数，默认值为 1.0，它表示采样间隔，也就是相邻两个样本点之间的时间或空间间隔
    frequency = numpy.fft.fftfreq(M, dt)
    # frequency = numpy.fft.fftfreq(N, dt)

    # 上述生成的是频率轴，我们希望得到角频率轴
    frequency_omega = 2.0 * math.pi * frequency
    
    # 幅度校正（补偿窗函数损耗）
    if apply_window:
        window_gain = numpy.mean(window)  # 窗函数增益
        frequency_spectrum /= window_gain
    
    ## return frequency_omega, frequency_spectrum
    return frequency, frequency_omega, frequency_spectrum

    ## 窗类型	特点	                          适用场景
    ## Hann	主瓣宽，旁瓣衰减快       黑洞准周期振荡(QPO)
    ## Tukey	可调平台区	     引力波啁啾信号
    ## Blackman	旁瓣抑制最优	     高动态范围数据

####################################################################################


####################################################################################

def frequency_filter_integration(omega, frequency_spectrum, omega0):

    '''
    ## 修改 omega < omega0 的部分
    omega_filter = numpy.where(omega < omega0, omega0, omega)
    '''

    ## 修改 |omega| < omega0 的部分
    
    # 根据符号生成替换值：正数→ omega0 ，负数→ - omega0
    replacements = numpy.where(omega >= 0, omega0, -omega0)
    
    # 应用替换
    omega_filter = numpy.where( numpy.abs(omega) < omega0, replacements, omega )


    ## 频域积分的被积函数
    ## 注意：时域的卷积相当于频域的乘积
    ## 
    frequency_integration = - frequency_spectrum / (omega_filter)**2

    return frequency_integration

####################################################################################


####################################################################################

## 该函数将 |omega| < omega0 的部分替换为 omega0

def omega_filter(omega, omega0):

    ## 修改 |omega| < omega0 的部分
    
    # 根据符号生成替换值：正数→ omega0 ，负数→ - omega0
    replacements = numpy.where(omega >= 0, omega0, -omega0)
    
    # 应用替换
    # omega_filter = numpy.where(mask, replacements, omega)
    omega_filter = numpy.where( numpy.abs(omega) < omega0, replacements, omega )

    return omega_filter
    
####################################################################################


####################################################################################

## 该函数用于进行逆傅里叶变换

## 需要的输入参数：
## omega : 频率轴（Hz）
## F_omega : 复数频域数据
## sampling_factor : 设置采样倍数
## original_zero_pad_factor : 原数据的零填充倍数
    
## 该函数的返回：
## t : 时间轴
## reconstructed : 重建的时域信号


def inverse_fourier_transform(omega, F_omega, sampling_factor=2, original_zero_pad_factor=4):
    
    # 计算采样参数
    N = len(F_omega)
    domega = omega[1] - omega[0]  # 频率分辨率

    # 设置采样率
    # 为了避免频谱混叠，采样率必须满足奈奎斯特采样定理，即采样率至少是信号最高频率的两倍
    if sampling_factor > 2: 
        sampling_rate_omega = sampling_factor * omega.max()
    else:
        sampling_rate_omega = 2.0 * omega.max()
    # 由于这里输入的是频率，因此我们乘以 2 pi
    ## dt = 1.0 / sampling_rate_omega
    frequency = omega / (2.0*math.pi)
    dt = 2 * math.pi / sampling_rate_omega  
    
    '''
    # 零频分量检查
    if not numpy.isclose(omega[0], 0):
        warnings.warn("频域数据未包含零频分量，可能导致直流偏移")
    '''
    
    # 执行逆傅里叶变换（与正变换选取相同的归一化模式）
    reconstructed_signal = numpy.fft.ifft(F_omega, norm='ortho') 
    # 注意，numpy.fft.ifft 本身就是归一化的，不需要乘以 N
    
    ## 生成时间轴
    ## 如果有零填充，则重建零填充数据
    if (original_zero_pad_factor > 1):
        N0 = N // original_zero_pad_factor
        t = numpy.arange(0, N0*dt, dt)
        reconstructed_signal2 = reconstructed_signal[:N0]
    ## 如果无零填充
    else:
        t = numpy.arange(0, N*dt, dt)
        reconstructed_signal2 = reconstructed_signal
    
    # 处理实数信号
    if numpy.allclose(numpy.imag(reconstructed_signal2), 0):
        reconstructed_signal3 = numpy.real(reconstructed_signal2)
    
    return t[:len(reconstructed_signal3)], reconstructed_signal3


####################################################################################


####################################################################################

# 该函数利用解析信号法 (Hilbert变换法)计算信号的瞬时频率

def instantaneous_frequency(signal, sampling_rate):
    """
    计算信号的瞬时频率
    :param signal: 输入的时间采样信号
    :param sampling_rate: 采样率
    :return: 时间序列和对应的瞬时频率序列
    """
    analytic_signal = scipy.signal.hilbert(signal)
    phase = numpy.unwrap(numpy.angle(analytic_signal))
    time = numpy.arange(len(signal)) / sampling_rate
    frequency = numpy.gradient(phase, time) / (2 * numpy.pi)
    return time, frequency

def get_frequency_at_t1(signal, sampling_rate, t1):
    """
    获取t1时刻的瞬时频率
    :param signal: 输入的时间采样信号
    :param sampling_rate: 采样率
    :param t1: 目标时刻
    :return: t1时刻的瞬时频率
    """
    time, freq = instantaneous_frequency(signal, sampling_rate)
    index = numpy.argmin(numpy.abs(time - t1))
    return freq[index]
    
    
####################################################################################



####################################################################################

## 该函数对引力波波形 h 画图

## 需要的输入参数：
## outdir             数据的文件夹地址
## figure_outdir      需要画图的文件夹地址
## detector_number_i  探测器序号
## total_mass         系统总质量

def generate_gravitational_wave_amplitude_plot( outdir, figure_outdir, detector_number_i ):


    print(                                                       )
    print( " 对引力波强度 h 进行画图 "                              )
    print( " 对第 ", detector_number_i, " 个探测器半径数据进行画图 " )
    print(                                                       )

    # 打开文件路径
    file0 = os.path.join(outdir, "bssn_psi4.dat")
    
    print( " 对应数据文件为 ", file0 )
    
    # 读取整个文件数据，假设数据是以空格分隔的浮点数
    data = numpy.loadtxt(file0)
    
    # 取出 phi4 文件中各列的数据
    time                 = data[:,0]
    psi4_l2m2m_real      = data[:,1]
    psi4_l2m2m_imaginary = data[:,2]
    psi4_l2m1m_real      = data[:,3]
    psi4_l2m1m_imaginary = data[:,4]
    psi4_l2m0_real       = data[:,5]
    psi4_l2m0_imaginary  = data[:,6]
    psi4_l2m1_real       = data[:,7]
    psi4_l2m1_imaginary  = data[:,8]
    psi4_l2m2_real       = data[:,9]
    psi4_l2m2_imaginary  = data[:,10]
    
    # 报错
    # 这里 file0 只是个文件名，不涉及 file.open 操作
    # file0.close()
    
    # python 中除法会返回浮点数，因此这里设置为整除
    length = len(time) // input_data.Detector_Number 
    
    time2                 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2m_real2      = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2m_imaginary2 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1m_real2      = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1m_imaginary2 = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m0_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m0_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m1_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2_real2       = numpy.zeros( (input_data.Detector_Number, length) )
    psi4_l2m2_imaginary2  = numpy.zeros( (input_data.Detector_Number, length) )
    
    # 将数据拆分为各探测器半径对应的数据
    for i in range(input_data.Detector_Number):
        for j in range(length):
            time2[i,j]                 = time[                 j*input_data.Detector_Number + i ]
            psi4_l2m2m_real2[i,j]      = psi4_l2m2m_real[      j*input_data.Detector_Number + i ]
            psi4_l2m2m_imaginary2[i,j] = psi4_l2m2m_imaginary[ j*input_data.Detector_Number + i ]
            psi4_l2m1m_real2[i,j]      = psi4_l2m1m_real[      j*input_data.Detector_Number + i ]
            psi4_l2m1m_imaginary2[i,j] = psi4_l2m1m_imaginary[ j*input_data.Detector_Number + i ]
            psi4_l2m0_real2[i,j]       = psi4_l2m0_real[       j*input_data.Detector_Number + i ]
            psi4_l2m0_imaginary2[i,j]  = psi4_l2m0_imaginary[  j*input_data.Detector_Number + i ]
            psi4_l2m1_real2[i,j]       = psi4_l2m1_real[       j*input_data.Detector_Number + i ]
            psi4_l2m1_imaginary2[i,j]  = psi4_l2m1_imaginary[  j*input_data.Detector_Number + i ]
            psi4_l2m2_real2[i,j]       = psi4_l2m2_real[       j*input_data.Detector_Number + i ]
            psi4_l2m2_imaginary2[i,j]  = psi4_l2m2_imaginary[  j*input_data.Detector_Number + i ]

    
    ## 对引力波数据 Psi4 进行离散傅里叶变换
    ## l=2 m=-2 频谱
    psi4_l2m2m_real_frequency, psi4_l2m2m_real_omega, psi4_l2m2m_real_omega_spectrem                                                             \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m2m_real2[detector_number_i],      apply_window=True, zero_pad_factor=4 )
    psi4_l2m2m_imaginary_frequency, psi4_l2m2m_imaginary_omega, psi4_l2m2m_imaginary_omega_spectrem                                              \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m2m_imaginary2[detector_number_i], apply_window=True, zero_pad_factor=4 )
    ## l=2 m=-1 频谱
    psi4_l2m1m_real_frequency, psi4_l2m1m_real_omega, psi4_l2m1m_real_omega_spectrem                                                             \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m1m_real2[detector_number_i],      apply_window=True, zero_pad_factor=4 )
    psi4_l2m1m_imaginary_frequency, psi4_l2m1m_imaginary_omega, psi4_l2m1m_imaginary_omega_spectrem                                              \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m1m_imaginary2[detector_number_i], apply_window=True, zero_pad_factor=4 )
    ## l=2 m=0 频谱
    psi4_l2m0_real_frequency, psi4_l2m0_real_omega, psi4_l2m0_real_omega_spectrem                                                                \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m0_real2[detector_number_i],       apply_window=True, zero_pad_factor=4 )
    psi4_l2m0_imaginary_frequency, psi4_l2m0_imaginary_omega, psi4_l2m0_imaginary_omega_spectrem                                                 \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m0_imaginary2[detector_number_i],  apply_window=True, zero_pad_factor=4 )
    ## l=2 m=1 频谱
    psi4_l2m1_real_frequency, psi4_l2m1_real_omega, psi4_l2m1_real_omega_spectrem                                                                \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m1_real2[detector_number_i],       apply_window=True, zero_pad_factor=4 )
    psi4_l2m1_imaginary_frequency, psi4_l2m1_imaginary_omega, psi4_l2m1_imaginary_omega_spectrem                                                 \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m1_imaginary2[detector_number_i],  apply_window=True, zero_pad_factor=4 )
    ## l=2 m=2 频谱
    psi4_l2m2_real_frequency, psi4_l2m2_real_omega, psi4_l2m2_real_omega_spectrem                                                                \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m2_real2[detector_number_i],       apply_window=True, zero_pad_factor=4 )
    psi4_l2m2_imaginary_frequency, psi4_l2m2_imaginary_omega, psi4_l2m2_imaginary_omega_spectrem                                                 \
        = compute_frequency_spectrum( time2[detector_number_i], psi4_l2m2_imaginary2[detector_number_i],  apply_window=True, zero_pad_factor=4 )


    # 根据输入数据推算出探测器距离
    Detector_Interval   = ( input_data.Detector_Rmax - input_data.Detector_Rmin ) / ( input_data.Detector_Number - 1 )
    Detector_Distance_R = input_data.Detector_Rmax - Detector_Interval * detector_number_i
    
    #################################################
    
    ## 计算系统总质量并输出 
    
    total_mass    = 0.0
    puncture_mass = numpy.zeros( input_data.puncture_number )
    
    ## 对于 Ansorg-TwoPuncture 的初值类型，对前两个黑洞质量归一化
    if ( input_data.Initial_Data_Method == "Ansorg-TwoPuncture" ):
        mass_ratio_Q = input_data.parameter_BH[0,0] / input_data.parameter_BH[1,0]
        BBH_M1 = mass_ratio_Q / ( 1.0 + mass_ratio_Q )
        BBH_M2 = 1.0          / ( 1.0 + mass_ratio_Q )
        for k in range( input_data.puncture_number ):
            if ( k == 0 ):
                puncture_mass[k] = BBH_M1 
            elif( k == 1 ):
                puncture_mass[k] = BBH_M2 
            else: 
                puncture_mass[k] = input_data.parameter_BH[k,0]
            total_mass += puncture_mass[k]
     
     ## 对于其它的初值类型，直接读入输入的 puncture 质量
    else:
        for k in range( input_data.puncture_number ):
            puncture_mass[k] = input_data.parameter_BH[k,0]
            total_mass += puncture_mass[k]
            
    ## 下面输出系统质量
    
    file1_path = os.path.join( figure_outdir, "tatal_mass.txt" )
    file1      = open( file1_path, "w" )
    
    for k in range( input_data.puncture_number ):
        print( f" mass[{k}] = {puncture_mass[k]} ", file=file1 )
        
    print( f" total mass = {total_mass} ", file=file1 )

    
    #################################################

    ## 设置乌龟坐标
    tortoise_R = Detector_Distance_R + 2.0 * total_mass * math.log( Detector_Distance_R / (2.0*total_mass) - 1.0)
    
    ## 设置初始时间
    t1 = tortoise_R

    ## 计算 Psi4 信号的瞬时频率
    ## instantaneous_frequency_psi4_l2m2_real = instantaneous_frequency( psi4_l2m2_real2[detector_number_i], len(psi4_l2m2_real2[detector_number_i]) )
    instantaneous_frequency_psi4_l2m2m_real = get_frequency_at_t1( psi4_l2m2m_real2[detector_number_i], len(psi4_l2m2m_real2[detector_number_i]), t1 )
    instantaneous_frequency_psi4_l2m1m_real = get_frequency_at_t1( psi4_l2m1m_real2[detector_number_i], len(psi4_l2m1m_real2[detector_number_i]), t1 )
    instantaneous_frequency_psi4_l2m0_real  = get_frequency_at_t1( psi4_l2m0_real2[detector_number_i],  len(psi4_l2m0_real2[detector_number_i]),  t1 )
    instantaneous_frequency_psi4_l2m1_real  = get_frequency_at_t1( psi4_l2m1_real2[detector_number_i],  len(psi4_l2m1_real2[detector_number_i]),  t1 )
    instantaneous_frequency_psi4_l2m2_real  = get_frequency_at_t1( psi4_l2m2_real2[detector_number_i],  len(psi4_l2m2_real2[detector_number_i]),  t1 )
    print( f" t - r* = 0 时刻瞬时频率 l=2 m=-2 phi4_real = {instantaneous_frequency_psi4_l2m2m_real:.2f} 1/M" )
    print( f" t - r* = 0 时刻瞬时频率 l=2 m=-1 phi4_real = {instantaneous_frequency_psi4_l2m1m_real:.2f} 1/M" )
    print( f" t - r* = 0 时刻瞬时频率 l=2 m=0  phi4_real = {instantaneous_frequency_psi4_l2m0_real:.2f}  1/M" )
    print( f" t - r* = 0 时刻瞬时频率 l=2 m=1  phi4_real = {instantaneous_frequency_psi4_l2m1_real:.2f}  1/M" )
    print( f" t - r* = 0 时刻瞬时频率 l=2 m=2  phi4_real = {instantaneous_frequency_psi4_l2m2_real:.2f}  1/M" )

    ## 根据瞬时频率添加频率截断条件
    frequency_cut_l2m2m = abs(instantaneous_frequency_psi4_l2m2m_real) * 1.2
    frequency_cut_l2m1m = abs(instantaneous_frequency_psi4_l2m1m_real) * 1.2
    frequency_cut_l2m0  = abs(instantaneous_frequency_psi4_l2m0_real)  * 1.2
    frequency_cut_l2m1  = abs(instantaneous_frequency_psi4_l2m1_real)  * 1.2
    frequency_cut_l2m2  = abs(instantaneous_frequency_psi4_l2m2_real)  * 1.2
    
    ## 添加频域滤波条件
    omega_cut_l2m2m = 2*math.pi / frequency_cut_l2m2m
    omega_cut_l2m1m = 2*math.pi / frequency_cut_l2m1m
    omega_cut_l2m0  = 2*math.pi / frequency_cut_l2m0
    omega_cut_l2m1  = 2*math.pi / frequency_cut_l2m1
    omega_cut_l2m2  = 2*math.pi / frequency_cut_l2m2
    
    ## 得到逆傅里叶变换的被积函数
    psi4_l2m2m_real_omega_integration      = frequency_filter_integration( psi4_l2m2m_real_omega,      psi4_l2m2m_real_omega_spectrem,      omega_cut_l2m2m )
    psi4_l2m2m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m2m_imaginary_omega, psi4_l2m2m_imaginary_omega_spectrem, omega_cut_l2m2m )
    psi4_l2m1m_real_omega_integration      = frequency_filter_integration( psi4_l2m1m_real_omega,      psi4_l2m1m_real_omega_spectrem,      omega_cut_l2m1m )
    psi4_l2m1m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m1m_imaginary_omega, psi4_l2m1m_imaginary_omega_spectrem, omega_cut_l2m1m )
    psi4_l2m0_real_omega_integration       = frequency_filter_integration( psi4_l2m0_real_omega,       psi4_l2m0_real_omega_spectrem,       omega_cut_l2m0  )
    psi4_l2m0_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m0_imaginary_omega,  psi4_l2m0_imaginary_omega_spectrem,  omega_cut_l2m0  )
    psi4_l2m1_real_omega_integration       = frequency_filter_integration( psi4_l2m1_real_omega,       psi4_l2m1_real_omega_spectrem,       omega_cut_l2m1  )
    psi4_l2m1_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m1_imaginary_omega,  psi4_l2m1_imaginary_omega_spectrem,  omega_cut_l2m1  )
    psi4_l2m2_real_omega_integration       = frequency_filter_integration( psi4_l2m2_real_omega,       psi4_l2m2_real_omega_spectrem,       omega_cut_l2m2  )
    psi4_l2m2_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m2_imaginary_omega,  psi4_l2m2_imaginary_omega_spectrem,  omega_cut_l2m2  )

    ## 在频域内进行逆傅里叶变换，得到引力波振幅
    ## l=2 m=-2 振幅
    time_grid_h_plus_l2m2m, GW_h_plus_l2m2m \
        = inverse_fourier_transform( psi4_l2m2m_real_omega, psi4_l2m2m_real_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  
    time_grid_h_cross_l2m2m, GW_h_cross_l2m2m \
        = inverse_fourier_transform( psi4_l2m2m_imaginary_omega, psi4_l2m2m_imaginary_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )
    ## l=2 m=-1 振幅
    time_grid_h_plus_l2m1m, GW_h_plus_l2m1m \
        = inverse_fourier_transform( psi4_l2m1m_real_omega, psi4_l2m1m_real_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  
    time_grid_h_cross_l2m1m, GW_h_cross_l2m1m \
        = inverse_fourier_transform( psi4_l2m1m_imaginary_omega, psi4_l2m1m_imaginary_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )
    ## l=2 m=0 振幅
    time_grid_h_plus_l2m0, GW_h_plus_l2m0 \
        = inverse_fourier_transform( psi4_l2m0_real_omega, psi4_l2m0_real_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  
    time_grid_h_cross_l2m0, GW_h_cross_l2m0 \
        = inverse_fourier_transform( psi4_l2m0_imaginary_omega, psi4_l2m0_imaginary_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )
    ## l=2 m=1 振幅
    time_grid_h_plus_l2m1, GW_h_plus_l2m1 \
        = inverse_fourier_transform( psi4_l2m1_real_omega, psi4_l2m1_real_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  
    time_grid_h_cross_l2m1, GW_h_cross_l2m1 \
        = inverse_fourier_transform( psi4_l2m1_imaginary_omega, psi4_l2m1_imaginary_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )
    ## l=2 m=2 振幅
    time_grid_h_plus_l2m2, GW_h_plus_l2m2 \
        = inverse_fourier_transform( psi4_l2m2_real_omega, psi4_l2m2_real_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  
    time_grid_h_cross_l2m2, GW_h_cross_l2m2 \
        = inverse_fourier_transform( psi4_l2m2_imaginary_omega, psi4_l2m2_imaginary_omega_integration, sampling_factor=2, original_zero_pad_factor=4 )  

    
    # 构造计算引力波振幅 h 的时间网格
    # time_max = max( time2[detector_number_i] )
    # time_grid = numpy.linspace( tortoise_R, time_max, 2000 )
    # time_grid_new = numpy.linspace( 0, time_max-tortoise_R, len(time_grid) )  ## 将时间减去乌龟坐标
    # l=2 m=-2
    time_grid_h_plus_l2m2m_new  = time_grid_h_plus_l2m2m - tortoise_R
    time_grid_h_cross_l2m2m_new = time_grid_h_cross_l2m2m - tortoise_R
    # l=2 m=-1
    time_grid_h_plus_l2m1m_new  = time_grid_h_plus_l2m1m - tortoise_R
    time_grid_h_cross_l2m1m_new = time_grid_h_cross_l2m1m - tortoise_R
    # l=2 m=0
    time_grid_h_plus_l2m0_new  = time_grid_h_plus_l2m0 - tortoise_R
    time_grid_h_cross_l2m0_new = time_grid_h_cross_l2m0 - tortoise_R
    # l=2 m=1
    time_grid_h_plus_l2m1_new  = time_grid_h_plus_l2m1 - tortoise_R
    time_grid_h_cross_l2m1_new = time_grid_h_cross_l2m1 - tortoise_R
    # l=2 m=2
    time_grid_h_plus_l2m2_new  = time_grid_h_plus_l2m2 - tortoise_R
    time_grid_h_cross_l2m2_new = time_grid_h_cross_l2m2 - tortoise_R

    plt.figure( figsize=(8,8) )                                   ## 这里 figsize 可以设定图形的大小
    plt.title( f" Gravitational Wave h   Detector Distence = { Detector_Distance_R } ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    plt.plot( time_grid_h_plus_l2m0_new, GW_h_plus_l2m0,  \
              color='red',    label="l=2 m=0 h+",                  linewidth=2 )
    plt.plot( time_grid_h_cross_l2m0_new, GW_h_cross_l2m0, \
              color='orange', label="l=2 m=0 hx",  linestyle='--', linewidth=2 )
    plt.plot( time_grid_h_plus_l2m1_new, GW_h_plus_l2m1,  \
              color='green',  label="l=2 m=1 h+",                  linewidth=2 )
    plt.plot( time_grid_h_cross_l2m1_new, GW_h_cross_l2m1, \
              color='cyan',   label="l=2 m=1 hx",  linestyle='--', linewidth=2 )
    plt.plot( time_grid_h_plus_l2m2_new, GW_h_plus_l2m2,  \
              color='black',  label="l=2 m=2 h+",                  linewidth=2 )
    plt.plot( time_grid_h_cross_l2m2_new, GW_h_cross_l2m2, \
              color='gray',   label="l=2 m=2 hx",  linestyle='--', linewidth=2 )
    plt.xlabel( "T [M]",          fontsize=16     )
    plt.ylabel( r"R*h",           fontsize=16     )
    plt.xlim( 0.0, max(time_grid_h_plus_l2m0_new) )
    plt.legend( loc='upper right'                 )
    plt.savefig( os.path.join(figure_outdir, "Gravitational_Wave_h_Detector_" + str(detector_number_i) + ".pdf") )
    
    print(                                                     )
    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " 对引力波强度 h 的画图完成 "                           )
    print(                                                     )
    
    '''
    # 以下为直接对Psi4积分，由于精度不够，已弃用
    # h = int_{0}^{t} dt' int_{0}^{t"} Psi4(t") dt"

    # 将各探测器数据进行插值，得到光滑函数
    # 这里使用三次样条插值
    psi4_l2m2m_real2_interpolation      = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m2m_real2[detector_number_i],      kind='cubic' )
    psi4_l2m2m_imaginary2_interpolation = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m2m_imaginary2[detector_number_i], kind='cubic' )
    psi4_l2m1m_real2_interpolation      = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m1m_real2[detector_number_i],      kind='cubic' )
    psi4_l2m1m_imaginary2_interpolation = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m1m_imaginary2[detector_number_i], kind='cubic' )
    psi4_l2m0_real2_interpolation       = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m0_real2[detector_number_i],       kind='cubic' )
    psi4_l2m0_imaginary2_interpolation  = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m0_imaginary2[detector_number_i],  kind='cubic' )
    psi4_l2m1_real2_interpolation       = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m1_real2[detector_number_i],       kind='cubic' )
    psi4_l2m1_imaginary2_interpolation  = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m1_imaginary2[detector_number_i],  kind='cubic' )
    psi4_l2m2_real2_interpolation       = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m2_real2[detector_number_i],       kind='cubic' )
    psi4_l2m2_imaginary2_interpolation  = scipy.interpolate.interp1d( time2[detector_number_i], psi4_l2m2_imaginary2[detector_number_i],  kind='cubic' )

    # 根据输入数据推算出探测器距离
    Detector_Interval   = ( input_data.Detector_Rmax - input_data.Detector_Rmin ) / ( input_data.Detector_Number - 1 )
    Detector_Distance_R = input_data.Detector_Rmax - Detector_Interval * detector_number_i

    # 设置乌龟坐标
    tortoise_R = Detector_Distance_R + 2.0 * total_mass * math.log( Detector_Distance_R / (2.0*total_mass) - 1.0)
    
    # 构造计算引力波振幅 h 的时间网格
    time_max = max( time2[detector_number_i] )
    time_grid = numpy.linspace( tortoise_R, time_max, 2000 )
    time_grid_new = numpy.linspace( 0, time_max-tortoise_R, len(time_grid) )  ## 将时间减去乌龟坐标

    GW_h_plus_l2m2m  = numpy.zeros( len(time_grid) )
    GW_h_cross_l2m2m = numpy.zeros( len(time_grid) )
    GW_h_plus_l2m1m  = numpy.zeros( len(time_grid) )
    GW_h_cross_l2m1m = numpy.zeros( len(time_grid) )
    GW_h_plus_l2m0   = numpy.zeros( len(time_grid) )
    GW_h_cross_l2m0  = numpy.zeros( len(time_grid) )
    GW_h_plus_l2m1   = numpy.zeros( len(time_grid) )
    GW_h_cross_l2m1  = numpy.zeros( len(time_grid) )
    GW_h_plus_l2m2   = numpy.zeros( len(time_grid) )
    GW_h_cross_l2m2  = numpy.zeros( len(time_grid) )

    # 通过数值积分求解 h = int_{0}^{t} dt' int_{0}^{t"} Psi4(t") dt" 
    # 该积分可以交换次序，化简为 h = int_{0}^{t} (t-t") Psi4(t") dt"
    def GW_h_plus_l2m2m_integrand(t, tmax):
        return psi4_l2m2m_real2_interpolation(t) * (tmax-t)
    def GW_h_cross_l2m2m_integrand(t, tmax):
        return psi4_l2m2m_imaginary2_interpolation(t) * (tmax-t)
    def GW_h_plus_l2m1m_integrand(t, tmax):
        return psi4_l2m1m_real2_interpolation(t) * (tmax-t)
    def GW_h_cross_l2m1m_integrand(t, tmax):
        return psi4_l2m1m_imaginary2_interpolation(t) * (tmax-t)
    def GW_h_plus_l2m0_integrand(t, tmax):
        return psi4_l2m0_real2_interpolation(t) * (tmax-t)
    def GW_h_cross_l2m0_integrand(t, tmax):
        return psi4_l2m0_imaginary2_interpolation(t) * (tmax-t)
    def GW_h_plus_l2m1_integrand(t, tmax):
        return psi4_l2m1_real2_interpolation(t) * (tmax-t)
    def GW_h_cross_l2m1_integrand(t, tmax):
        return psi4_l2m1_imaginary2_interpolation(t) * (tmax-t)
    def GW_h_plus_l2m2_integrand(t, tmax):
        return psi4_l2m2_real2_interpolation(t) * (tmax-t)
    def GW_h_cross_l2m2_integrand(t, tmax):
        return psi4_l2m2_imaginary2_interpolation(t) * (tmax-t)
    
    # 计算引力波振幅 h+ 和 hx
    # 在积分中用 lambda 函数重定义被积函数，使之变为单变量函数

    for j in range( len(time_grid) ):

        print( " j = ", j )

        GW_h_plus_l2m2m_integrand2 = lambda t: GW_h_plus_l2m2m_integrand(t, time_grid[j])
        ## 注意这里 scipy.integrate.quad 返回的是元组，第一个是积分值，第二个是误差
        GW_h_plus_l2m2m[j], err0 = scipy.integrate.quad( GW_h_plus_l2m2m_integrand2, 0.0, time_grid[j], limit=600 )
                                                   # epsabs=1e-8,  # 绝对误差
                                                   # limit=600 )    # 增加分段数 )

        GW_h_cross_l2m2m_integrand2 = lambda t: GW_h_cross_l2m2m_integrand(t, time_grid[j])
        GW_h_cross_l2m2m[j], err0 = scipy.integrate.quad( GW_h_cross_l2m2m_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_plus_l2m1m_integrand2 = lambda t: GW_h_plus_l2m1m_integrand(t, time_grid[j])
        GW_h_plus_l2m1m[j], err0 = scipy.integrate.quad( GW_h_plus_l2m1m_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_cross_l2m1m_integrand2 = lambda t: GW_h_cross_l2m1m_integrand(t, time_grid[j])
        GW_h_cross_l2m1m[j], err0 = scipy.integrate.quad( GW_h_cross_l2m1m_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_plus_l2m0_integrand2 = lambda t: GW_h_plus_l2m0_integrand(t, time_grid[j])
        GW_h_plus_l2m0[j], err0 = scipy.integrate.quad( GW_h_plus_l2m0_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_cross_l2m0_integrand2 = lambda t: GW_h_cross_l2m0_integrand(t, time_grid[j])
        GW_h_cross_l2m0[j], err0 = scipy.integrate.quad( GW_h_cross_l2m0_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_plus_l2m1_integrand2 = lambda t: GW_h_plus_l2m1_integrand(t, time_grid[j])
        GW_h_plus_l2m1[j], err0 = scipy.integrate.quad( GW_h_plus_l2m1_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_cross_l2m1_integrand2 = lambda t: GW_h_cross_l2m1_integrand(t, time_grid[j])
        GW_h_cross_l2m1[j], err0 = scipy.integrate.quad( GW_h_cross_l2m1_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_plus_l2m2_integrand2 = lambda t: GW_h_plus_l2m2_integrand(t, time_grid[j])
        GW_h_plus_l2m2[j], err0 = scipy.integrate.quad( GW_h_plus_l2m2_integrand2, 0.0, time_grid[j], limit=600 )

        GW_h_cross_l2m2_integrand2 = lambda t: GW_h_cross_l2m2_integrand(t, time_grid[j])
        GW_h_cross_l2m2[j], err0 = scipy.integrate.quad( GW_h_cross_l2m2_integrand2, 0.0, time_grid[j], limit=600 )
            
    # 对引力波振幅 h+ 和 hx 的计算完成

    # 下面进行画图
    plt.figure( figsize=(8,8) )                                   ## 这里 figsize 可以设定图形的大小
    plt.title( f" Gravitational Wave h   Detector Distence = { Detector_Distance_R } ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    plt.plot( time_grid_new, GW_h_plus_l2m0,  \
              color='red',    label="l=2 m=0 h+",                  linewidth=2 )
    plt.plot( time_grid_new, GW_h_cross_l2m0, \
              color='orange', label="l=2 m=0 hx",  linestyle='--', linewidth=2 )
    plt.plot( time_grid_new, GW_h_plus_l2m1,  \
              color='green',  label="l=2 m=1 h+",                  linewidth=2 )
    plt.plot( time_grid_new, GW_h_cross_l2m1, \
              color='cyan',   label="l=2 m=1 hx",  linestyle='--', linewidth=2 )
    plt.plot( time_grid_new, GW_h_plus_l2m2,  \
              color='black',  label="l=2 m=2 h+",                  linewidth=2 )
    plt.plot( time_grid_new, GW_h_cross_l2m2, \
              color='gray',   label="l=2 m=2 hx",  linestyle='--', linewidth=2 )
    plt.xlabel( "T [M]",          fontsize=16 )
    plt.ylabel( r"R*h",           fontsize=16 )
    plt.legend( loc='upper right'             )
    plt.savefig( os.path.join(figure_outdir, "Gravitational_Wave_h_Detector_" + str(detector_number_i) + ".pdf") )
    
    print(                                                     )
    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " 对引力波强度 h 的画图完成 "                           )
    print(                                                     )
    '''

    return

####################################################################################



####################################################################################

## 单独使用的例子

## outdir = "GW150914"
## generate_gravitational_wave_amplitude_plot(outdir, outdir, 0)

####################################################################################



