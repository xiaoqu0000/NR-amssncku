
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
import derivative_xiaoqu

# plt.rcParams['text.usetex'] = True  ## 在绘图中允许使用 latex 字体

#from scipy import signal
#import matplotlib.pyplot as plt

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

## 该函数用于进行逆傅里叶变换

## 需要的输入参数：
## omega : 频率轴（Hz）
## F_omega : 复数频域数据
## sampling_factor : 设置采样倍数
## original_length : 原始信号长度（可选）
    
## 该函数的返回：
## t : 时间轴
## reconstructed : 重建的时域信号


def inverse_fourier_transform(omega, F_omega, sampling_factor=2, original_zero_pad_factor=4):
## def inverse_fourier_transform(omega, F_omega, sampling_rate=None, original_length=None):
    
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

## 计算初始轨道的瞬时周期和瞬时频率

def estimate_orbit_frequency_at_t0( data, method="7-points 6-orders" ):

    ## 设定最大频率和最小频率
    frequency_min = 1.0
    frequency_max = 1000.0

    # 提取时间和坐标
    t = data[:,0]
    x = data[:,1]
    y = data[:,2]

    # 计算相位
    phase = numpy.arctan2(y, x)

    ## 如果数组长度足够，利用相位的变化率求出对应周期和频率
    if (len(t) >= 10 ):
        dphi_dt = derivative_xiaoqu.first_order_derivative_at_t0( phase, t, 5, method )
        omega     = abs(dphi_dt)
        frequency = omega / ( 2.0 * math.pi )
        period_T  = 1.0 / frequency
        print( " omega = ", omega )
        print( " period = ", period_T )
        if (frequency < frequency_min):
            frequency = frequency_min
            period_T  = 1.0 / frequency
            omega     = frequency_min * 2.0 * math.pi
        if (frequency > frequency_min):
            frequency = frequency_max
            period_T  = 1.0 / frequency
            omega     = frequency_max * 2.0 * math.pi
    ## 如果数组长度不够，利用最大频率赋值
    else:
        frequency = frequency_max
        period_T  = 1.0 / frequency
        omega     = frequency_max * 2.0 * math.pi

    return frequency, omega, period_T


####################################################################################

## 该函数对引力波波形 h 画图

def generate_gravitational_amplitude_plot( outdir, figure_outdir, detector_number_i, total_mass, frequency_cut ):


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

    # 设置乌龟坐标
    tortoise_R = Detector_Distance_R + 2.0 * total_mass * math.log( Detector_Distance_R / (2.0*total_mass) - 1.0)
    
    ## 添加频域滤波条件 
    t1 = tortoise_R

    ## instantaneous_frequency_psi4_l2m2_real = instantaneous_frequency( psi4_l2m2_real2[detector_number_i], len(psi4_l2m2_real2[detector_number_i]) )
    instantaneous_frequency_psi4_l2m2_real = get_frequency_at_t1( psi4_l2m2_real2[detector_number_i], len(psi4_l2m2_real2[detector_number_i]), t1)
    instantaneous_frequency_psi4_l2m1_real = get_frequency_at_t1( psi4_l2m1_real2[detector_number_i], len(psi4_l2m1_real2[detector_number_i]), t1)
    instantaneous_frequency_psi4_l2m0_real = get_frequency_at_t1( psi4_l2m0_real2[detector_number_i], len(psi4_l2m0_real2[detector_number_i]), t1)
    print(f" t - r* = 0 时刻瞬时频率 l=2 m=2 phi4_real = {instantaneous_frequency_psi4_l2m2_real:.2f} 1/M")
    print(f" t - r* = 0 时刻瞬时频率 l=2 m=1 phi4_real = {instantaneous_frequency_psi4_l2m1_real:.2f} 1/M")
    print(f" t - r* = 0 时刻瞬时频率 l=2 m=0 phi4_real = {instantaneous_frequency_psi4_l2m0_real:.2f} 1/M")

    
    ## 打开黑洞轨迹文件，估算初始时刻的轨道瞬时频率和周期
    file_BH = os.path.join(outdir, "bssn_BH.dat")
    
    ## 读取整个文件数据，假设数据是以空格分隔的浮点数
    BH_data = numpy.loadtxt(file_BH)

    ## 求出黑洞轨道瞬时频率和周期
    frequency_t0, omega_t0, period_T_t0 = estimate_orbit_frequency_at_t0( BH_data, method="5-points 4-orders" )
    frequency_t0, omega_t0, period_T_t0 = estimate_orbit_frequency_at_t0( BH_data, method="7-points 6-orders" )
    print(f" t - r* = 0 时刻轨道瞬时频率 = {frequency_t0:.2f} 1/M")
    print(f" t - r* = 0 时刻轨道瞬时周期 = {period_T_t0:.2f}    M")
    
    ## 添加频率截断条件
    period_cut = period_T_t0 
    omega_cut  = 2.0 * math.pi / period_cut
    period_cut = 2.0 * math.pi * ( (11*11*11)**0.5 ) 
    print( "period cut = ", period_cut )
    omega_cut  = ( 2.0 * math.pi / period_cut ) * 2
    
    ## 得到逆傅里叶变换的被积函数
    psi4_l2m2m_real_omega_integration      = frequency_filter_integration( psi4_l2m2m_real_omega,      psi4_l2m2m_real_omega_spectrem,      omega_cut )
    psi4_l2m2m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m2m_imaginary_omega, psi4_l2m2m_imaginary_omega_spectrem, omega_cut )
    psi4_l2m1m_real_omega_integration      = frequency_filter_integration( psi4_l2m1m_real_omega,      psi4_l2m1m_real_omega_spectrem,      omega_cut )
    psi4_l2m1m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m1m_imaginary_omega, psi4_l2m1m_imaginary_omega_spectrem, omega_cut )
    psi4_l2m0_real_omega_integration       = frequency_filter_integration( psi4_l2m0_real_omega,       psi4_l2m0_real_omega_spectrem,       omega_cut )
    psi4_l2m0_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m0_imaginary_omega,  psi4_l2m0_imaginary_omega_spectrem,  omega_cut )
    psi4_l2m1_real_omega_integration       = frequency_filter_integration( psi4_l2m1_real_omega,       psi4_l2m1_real_omega_spectrem,       omega_cut )
    psi4_l2m1_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m1_imaginary_omega,  psi4_l2m1_imaginary_omega_spectrem,  omega_cut )
    psi4_l2m2_real_omega_integration       = frequency_filter_integration( psi4_l2m2_real_omega,       psi4_l2m2_real_omega_spectrem,       omega_cut )
    psi4_l2m2_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m2_imaginary_omega,  psi4_l2m2_imaginary_omega_spectrem,  omega_cut )

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
    
    print(                                                       )
    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " 对引力波强度 h 的画图完成 "                           )
    print(                                                       )
    
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
    
    print(                                                       )
    print( " 第 ", detector_number_i, " 个探测器半径数据画图完成 " )
    print( " 对引力波强度 h 的画图完成 "                           )
    print(                                                       )
    '''

    return

####################################################################################




####################################################################################

## 这是一个傅里叶变换的测试

def test_fourier_transform( outdir, figure_outdir, detector_number_i, total_mass, frequency_cut ):


    ##############################################################################

    ## 对给定信号的测试

    ##############################################################################


    file1_path = os.path.join( outdir, "signal_text.txt"  )
    file2_path = os.path.join( outdir, "signal_text2.txt" )
    file1 = open( file1_path, "w" )
    file2 = open( file2_path, "w" )

    # 示例数据
    sampling_rate = 2001

    t_signal = numpy.linspace(0.0, 1.0, sampling_rate)
    signal   = numpy.sin( 2 * math.pi * 50 * t_signal )  # 50Hz的正弦信号
    dt       = t_signal[1] - t_signal[0]
    fs       = 1.0 / dt
    
    ## 计算频率轴
    frequency_fft = numpy.fft.fftfreq(sampling_rate, dt)
    omega_fft     = 2.0 * math.pi * frequency_fft
    
    t1 = 0.2
    freq_at_t1 = get_frequency_at_t1(signal, sampling_rate, t1)
    print(f" t = {t1} 时刻的瞬时频率: {freq_at_t1} Hz ")

    signal_omega = numpy.fft.fft(  signal,       norm='ortho' )     # 正变换
    signal_recon = numpy.fft.ifft( signal_omega, norm='ortho' )     # 逆变换

    frequency_xiaoqu, omega_xiaoqu, signal_omega_xiaoqu = compute_frequency_spectrum(t_signal, signal, apply_window=False, zero_pad_factor=4)

    print(                                      file=file1 )
    print( " omega    = ", omega_xiaoqu,        file=file1 )
    print(                                      file=file1 )
    print( " f(omega) = ", signal_omega_xiaoqu, file=file1 )
    print(                                      file=file1 )

    print( " i  omega_i   f(omega_i)  omega_xiaoqu_i  f(omega_xiaoqu_i) ", file=file2  )
    for i in range( len(omega_fft) ):
        print( format(i," 3d"),               ""*2,                                                                                                           \
               format(omega_fft[i],".5e"),    ""*2, format(signal_omega[i].imag,".5e"),        "+ (", format(signal_omega[i].imag,".5e"),        ")j", ""*2,  \
               format(omega_xiaoqu[i],".5e"), ""*2, format(signal_omega_xiaoqu[i].real,".5e"), "+ (", format(signal_omega_xiaoqu[i].imag,".5e"), ")j", ""*2,  \
               file=file2 )

    time_xiaoqu, signal_recon_xiaoqu = inverse_fourier_transform(omega_xiaoqu, signal_omega_xiaoqu, sampling_factor=4, original_zero_pad_factor=4)
    signal_recon_xiaoqu2 = signal_recon_xiaoqu[:sampling_rate]  ## 去掉零填充的信号，只保留有效信号
    
    print(                                           file=file1 )
    print( " signal        = ", signal,              file=file1 )
    print(                                           file=file1 )
    print( " signal_xiaoqu = ", signal_recon_xiaoqu, file=file1 )
    print(                                           file=file1 )

    ## 取出前一半数据（正频部分）进行画图
    ## 去掉虚数部分
    half_length              = len(omega_fft) // 2
    half_length2             = len(omega_xiaoqu) // 2
    omega_fft_half           = numpy.real(omega_fft)[:half_length]
    signal_omega_half        = numpy.real(signal_omega)[:half_length]
    omega_xiaoqu_half        = numpy.real(omega_xiaoqu)[:half_length2]
    signal_omega_xiaoqu_half = numpy.real(signal_omega_xiaoqu)[:half_length2]

    plt.figure( figsize=(8,8) )   ## 这里 figsize 可以设定图形的大小
    plt.plot( numpy.real(omega_fft),    numpy.real(signal_omega),        color='red',   label=r"$\omega(X)$",                       linewidth=2 )
    plt.plot( numpy.real(omega_xiaoqu), numpy.real(signal_omega_xiaoqu), color='green', label=r"$\omega_{xiaoqu}$", linestyle='--', linewidth=2 )
    plt.xlabel( "X",         fontsize=16  )
    plt.ylabel( r"$\omega$", fontsize=16  )
    plt.xlim(   -500.0, 500.0     )
    plt.legend( loc='upper right' )
    plt.savefig( os.path.join(figure_outdir, "test0a.pdf") )
    plt.close()

    plt.figure( figsize=(8,8) )   ## 这里 figsize 可以设定图形的大小
    ## plt.plot( numpy.real(omega_fft),    numpy.real(signal_omega),        color='red',   label=r"$\omega(X)$",                       linewidth=2 )
    ## plt.plot( numpy.real(omega_xiaoqu), numpy.real(signal_omega_xiaoqu), color='green', label=r"$\omega_{xiaoqu}$", linestyle='--', linewidth=2 )
    plt.plot( omega_fft_half,    signal_omega_half,        color='red',   label=r"$\omega(X)$",                       linewidth=2 )
    plt.plot( omega_xiaoqu_half, signal_omega_xiaoqu_half, color='green', label=r"$\omega_{xiaoqu}$", linestyle='--', linewidth=2 )
    plt.xlabel( "X",         fontsize=16  )
    plt.ylabel( r"$\omega$", fontsize=16  )
    plt.xlim(   0.0, 500.0        )
    plt.legend( loc='upper right' )
    plt.savefig( os.path.join(figure_outdir, "test0b.pdf") )
    plt.close()

    plt.figure( figsize=(8,8) )   ## 这里 figsize 可以设定图形的大小
    plt.plot( t_signal, signal,               color='red',    label="signal",                                    linewidth=2 )
    plt.plot( t_signal, signal_recon,         color='orange', label="signal_reconstruct_ifft",   linestyle='--', linewidth=2 )
    plt.plot( t_signal, signal_recon_xiaoqu2, color='green',  label="signal_reconstruct_xiaoqu", linestyle='--', linewidth=2 )
    plt.xlabel( "T", fontsize=16  )
    plt.ylabel( "f", fontsize=16  )
    plt.legend( loc='upper right' )
    plt.savefig( os.path.join(figure_outdir, "test1.pdf") )
    plt.close()

    # 正确的归一化关系
    x_random   = numpy.random.rand( 501 )
    t_x_random = numpy.linspace(0.0, 1.0, 501, endpoint=False)
    # print( x[1] )

    x_random_omega = numpy.fft.fft(  x_random,       norm='ortho' )     # 正变换
    x_random_recon = numpy.fft.ifft( x_random_omega, norm='ortho' )     # 逆变换

    print(                                 file=file1 )
    print( "x_random = ",  x_random,       file=file1 )
    print(                                 file=file1 )
    print( "x_random2 = ", x_random_recon, file=file1 )
    print(                                 file=file1 )

    frequency_xiaoqu, omega_xiaoqu, x_random_omega_xiaoqu = compute_frequency_spectrum(t_x_random, x_random, apply_window=False, zero_pad_factor=4)
    print(                                              file=file1 )
    print( "x_random_omega = ",  x_random_omega,        file=file1 )
    print(                                              file=file1 )
    print( "x_random_omega2 = ", x_random_omega_xiaoqu, file=file1 )
    print(                                              file=file1 )

    t_xiaoqu, x_random_recon_xiaoqu = inverse_fourier_transform(omega_xiaoqu, x_random_omega_xiaoqu, sampling_factor=4, original_zero_pad_factor=4)
    x_random_recon_xiaoqu2 = x_random_recon_xiaoqu[:501]  ## 去掉零填充的信号，只保留有效信号

    print(                                              file=file1 )
    print( "x_random        = ", x_random,              file=file1 )
    print(                                              file=file1 )
    print( "x_random_xiaoqu = ", x_random_recon_xiaoqu, file=file1 )
    print(                                              file=file1 )

    file1.close()

    plt.figure( figsize=(8,8) )   ## 这里 figsize 可以设定图形的大小
    plt.plot( t_x_random, x_random,               color='red',    label="x_random",                                    linewidth=2 )
    plt.plot( t_x_random, x_random_recon,         color='orange', label="x_random_reconstruct_ifft",   linestyle='--', linewidth=2 )
    plt.plot( t_x_random, x_random_recon_xiaoqu2, color='green',  label="x_random_reconstruct_xiaoqu", linestyle='--', linewidth=2 )
    plt.xlabel( "T", fontsize=16  )
    plt.ylabel( "X", fontsize=16  )
    plt.legend( loc='upper right' )
    plt.savefig( os.path.join(figure_outdir, "test2.pdf") )
    plt.close()
    
    
    ##############################################################################

    ## 对引力波 Psi4 信号的测试

    ##############################################################################

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
    
    
    '''
    psi4_l2m2m_real2_omega          = numpy.empty(input_data.Detector_Number) 
    psi4_l2m2m_real2_omega_filter   = numpy.empty(input_data.Detector_Number)
    psi4_l2m2m_real2_omega_spectrem = numpy.empty(input_data.Detector_Number)
    psi4_l2m2m_real2_omega_integration = numpy.empty(input_data.Detector_Number)
    '''
    
    ## 测试代码
    # 正确的归一化关系
    '''
    x = psi4_l2m2_real2[detector_number_i]
    print( "psi4=", x )

    X = numpy.fft.fft(x, norm='ortho')            # 正变换
    print( "psi4(omega) = ", X)

    x_recon = numpy.fft.ifft(X, norm='ortho')     # 逆变换
    print( "x  = ", x       )
    print( "x2 = ", x_recon )
    '''
    
    ## 添加频率截断条件
    omega0 = 2*math.pi / frequency_cut
    
    ## 得到逆傅里叶变换的被积函数
    psi4_l2m2m_real_omega_integration      = frequency_filter_integration( psi4_l2m2m_real_omega,      psi4_l2m2m_real_omega_spectrem,      omega0 )
    psi4_l2m2m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m2m_imaginary_omega, psi4_l2m2m_imaginary_omega_spectrem, omega0 )
    psi4_l2m1m_real_omega_integration      = frequency_filter_integration( psi4_l2m1m_real_omega,      psi4_l2m1m_real_omega_spectrem,      omega0 )
    psi4_l2m1m_imaginary_omega_integration = frequency_filter_integration( psi4_l2m1m_imaginary_omega, psi4_l2m1m_imaginary_omega_spectrem, omega0 )
    psi4_l2m0_real_omega_integration       = frequency_filter_integration( psi4_l2m0_real_omega,       psi4_l2m0_real_omega_spectrem,       omega0 )
    psi4_l2m0_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m0_imaginary_omega,  psi4_l2m0_imaginary_omega_spectrem,  omega0 )
    psi4_l2m1_real_omega_integration       = frequency_filter_integration( psi4_l2m1_real_omega,       psi4_l2m1_real_omega_spectrem,       omega0 )
    psi4_l2m1_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m1_imaginary_omega,  psi4_l2m1_imaginary_omega_spectrem,  omega0 )
    psi4_l2m2_real_omega_integration       = frequency_filter_integration( psi4_l2m2_real_omega,       psi4_l2m2_real_omega_spectrem,       omega0 )
    psi4_l2m2_imaginary_omega_integration  = frequency_filter_integration( psi4_l2m2_imaginary_omega,  psi4_l2m2_imaginary_omega_spectrem,  omega0 )
    
    ## 测试代码
    ## 画傅里叶频谱
    psi4_l2m2_real_omega_filter = omega_filter(psi4_l2m2_real_omega, omega0)
    '''
    print( psi4_l2m2_real_omega ) 
    print( psi4_l2m2_real_omega_spectrem )
    print( psi4_l2m2_real_omega_filter ) 
    print( numpy.real(psi4_l2m2_real_omega_spectrem) )
    '''

    plt.figure( figsize=(8,8) )                                   ## 这里 figsize 可以设定图形的大小
    plt.title( f" Gravitational Wave h  ", fontsize=18 )   ## 这里 fontsize 可以设定文字大小
    ## plt.plot( time_grid_new, GW_h_plus_l2m0,  \
    ##           color='red',    label="l=2 m=0 h+",                  linewidth=2 )
    ## plt.plot( time_grid_new, GW_h_cross_l2m0, \
    ##           color='orange', label="l=2 m=0 hx",  linestyle='--', linewidth=2 )
    ## plt.plot( time_grid_new, GW_h_plus_l2m1,  \
    ##           color='green',  label="l=2 m=1 h+",                  linewidth=2 )
    ## plt.plot( time_grid_new, GW_h_cross_l2m1, \
    ##           color='cyan',   label="l=2 m=1 hx",  linestyle='--', linewidth=2 )
    plt.plot( numpy.log(psi4_l2m2_real_omega), numpy.log(numpy.real(psi4_l2m2_real_omega_spectrem)), \
              color='black',  label="l=2 m=2 $\psi$($\omega$)",                       linewidth=2 )
    plt.plot( numpy.log(psi4_l2m2_real_omega), numpy.log(psi4_l2m2_real_omega_integration),          \
              color='blue',  label=r"l=2 m=2 -$\psi$($\omega$)/$\tilde{\omega}^{2}$", linewidth=2 )
    ## plt.plot( psi4_l2m2_real_omega, numpy.real(psi4_l2m2_real_omega_spectrem),  \
    ##          color='black',  label="l=2 m=2 psi(omega)",                  linewidth=2 )
    ## plt.plot( time_grid_new, GW_h_cross_l2m2, \
    ##           color='gray',   label="l=2 m=2 hx",  linestyle='--', linewidth=2 )
    plt.xlabel( r"log($\omega$) [$M^{-1}$]",   fontsize=16 )
    plt.ylabel( r"log($\psi$($\omega$))", fontsize=16 )
    plt.legend( loc='upper right'                     )
    plt.savefig( os.path.join(figure_outdir, "Psi_omega_" + str(detector_number_i) + ".pdf") )
    
    return


####################################################################################


####################################################################################

# 单独使用的例子

outdir = "BBH_q=1"
# outdir = "BBH_q=1_new"
# generate_puncture_orbit_plot(outdir, outdir)
generate_gravitational_amplitude_plot(outdir, outdir, 10, 1.0, 120)
test_fourier_transform(outdir, outdir, 10, 1.0, 120)

####################################################################################


