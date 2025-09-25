# AMSS-NCKU 介绍

#### AMSS-NCKU 能做些什么

AMSS-NCKU是一个中国开发的数值相对论程序，用于对爱因斯坦方程进行数值求解，计算引力场随时间的变化。
AMSS-NCKU使用有限差分方法，通过自适应网格细化技术，来实现爱因斯坦方程的数值求解。
目前，AMSS-NCKU能够成功地处理双黑洞系统、多个黑洞系统，计算这些系统的时间演化，求解这些过程中释放的引力波

#### AMSS-NCKU 的发展历程

2008年，AMSS-NCKU程序开发成功，实现了BSSN方程的数值求解，用于双黑洞和多个黑洞系统。

2013年，AMSS-NCKU实现了Z4C方程的数值求解，提高了计算的精确度。

2015年，针对BSSN方程的求解，AMSS-NCKU实现了CPU和GPU混合计算，获得了7倍硬件加速。

2024年，我们为AMSS-NCKU开发了Python操作接口，方便用户使用和程序后续开发

#### AMSS-NCKU 的开发团队

曹周键（北京师范大学、中国科学院数学与系统科学研究院、中国科学院大学杭州高等科学研究院）

游辉樟（成功大学）

刘润球（中国科学院数学与系统科学研究院）

都志辉（清华大学）

季力伟（罗彻斯特理工学院）

赵志超（中国农业大学）

乔琛凯（重庆理工大学）

余瑞彬（曾在读学生）

林俊玉（曾在读学生）


#### 运行 AMSS-NCKU 所需要的依赖包

这里以 Ubuntu22.04 系统为例进行介绍

1.  安装需要的 C++/Fortran/Cuda 编译器

    sudo apt-get install gcc

    sudo apt-get install gfortran

    sudo apt-get install make

    sudo apt-get install build-essential

    sudo apt-get install nvidia-cuda-toolkit

2.  安装 MPI 工具

    sudo apt install openmpi-bin

    sudo apt install libopenmpi-dev

3.  安装 Python

    sudo apt-get install python3

    sudo apt-get install python3-pip

4.  安装所需要的 Python 库

    pip install numpy

    pip install scipy

    pip install matplotlib

    pip install SymPy

    pip install opencv-python-full

    pip install notebook

    pip install torch

5.  安装OpenCV工具

    sudo apt-get install libopencv-dev

    sudo apt-get install python-opencv

#### 如何使用 AMSS-NCKU

0.  修改编译相关的设置
    修改 AMSS_NCKU_source 目录下的 makefile.inc 文件，根据所用计算机来进行设置
    Ununtu22.04系统已经设置完成，不用修改

1.  进入AMSS-NCKU代码文件夹，修改输入
    输入设置在AMSS_NCKU_Input.py文件中，修改该文件并保存

2.  启动AMSS-NCKU
    控制台运行 python AMSS_NCKU_Program.py 或 python3 AMSS_NCKU_Program.py）

#### 更新记录

....

#### 重要提示

由于测试有限，程序难免会有小bug存在

进行一个实际的演化运行的时间比较久，为了避免中运行中途出现bug（例如运行结束后的绘图），可先将输入文件中的最终演化时间调为5，进行测试。
如果能成功运行并无报错，再调整输入文件到实际的最终演化时间（大约1000多），这样可以减少不必要的计算资源浪费。

请根据自身计算机设置计算资源（输入文件中设定MPI进程数目）

当前版本下双星准圆轨道轨道后牛顿计算的轨道动量初始值误差较大，导致自旋角动量较大的双星系统有时候不够圆（无自旋的情况稍好一些）。用户可在输入文件中选择Manually来自己设置更精确的双星动量初值。
稍后一段时间会发布新的版本提高准圆轨道初始动量的准确度。
