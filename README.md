# AMSS-NCKU 

#### AMSS-NCKU 能做些什么

AMSS-NCKU是一个中国开发的数值相对论程序，用于对爱因斯坦方程进行数值求解，计算引力场随时间的变化。

AMSS-NCKU使用有限差分方法，通过自适应网格细化技术，来实现爱因斯坦方程的数值求解。

目前，AMSS-NCKU能够成功地处理双黑洞系统、多个黑洞系统，计算这些系统的时间演化，求解这些过程中释放的引力波。

#### What can AMSS-NCKU do

AMSS - NCKU is a numerical relativity program developed in China, which is used to numerically solve Einstein's equations and calculate the change of the gravitational field over time. 

AMSS - NCKU uses the finite difference method and the adaptive mesh refinement technique to achieve the numerical solution of Einstein's equations. 

Currently, AMSS - NCKU can successfully handle binary black hole systems and multiple black hole systems, calculate the time evolution of these systems, and solve the gravitational waves released during these processes.

#### AMSS-NCKU 的发展历程

2008年，AMSS-NCKU程序开发成功，实现了BSSN方程的数值求解，用于双黑洞和多个黑洞系统。

2013年，AMSS-NCKU实现了Z4C方程的数值求解，提高了计算的精确度。

2015年，针对BSSN方程的求解，AMSS-NCKU实现了CPU和GPU混合计算。

2024年，我们为AMSS-NCKU开发了Python操作接口，方便用户使用和程序后续开发

#### The Development of AMSS-NCKU

In 2008, the AMSS-NCKU code was successfully developed, enabling the numerical simulation for binary black hole and multiple black hole systems via the BSSN equations .

In 2013, AMSS-NCKU achieved the numerical simulation for black hole systems via the Z4C equations, greatly improving the accuracy of calculation.

In 2015, AMSS-NCKU implemented hybrid CPU and GPU computing for the BSSN equations, improving the computational efficiency.

In 2024, we developed a Python operation interface for AMSS-NCKU to facilitate the freshman users and subsequent development.

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

左毅（学生）

#### Authors of AMSS-NCKU

Cao, Zhoujian (Beijing Normal University; Academy of Mathematics and Systems Science, Chinese Academy of Sciences; Hangzhou Institute for Advanced Study, University of Chinese Academy of Sciences)

Yo, Hwei-Jang (National Cheng Kung University)

Liu, Runqiu (Academy of Mathematics and Systems Science, Chinese Academy of Sciences)

Du, Zhihui (Tsinghua University)

Ji, Liwei (Rochester Institute of Technology)

Zhao, Zhichao (China Agricultural University)

Qiao, Chenkai (Chongqing University of Technology)

Yu, Jui-Ping (Former student)

Lin, Chun-Yu (Former student)

Zuo, Yi (Student)


#### 运行 AMSS-NCKU 所需要的依赖包

这里以 Ubuntu22.04 系统为例进行介绍

1.  安装需要的 C++/Fortran/Cuda 编译器

    $ sudo apt-get install gcc

    $ sudo apt-get install gfortran

    $ sudo apt-get install make

    $ sudo apt-get install build-essential

    $ sudo apt-get install nvidia-cuda-toolkit

2.  安装 MPI 工具

    $ sudo apt install openmpi-bin

    $ sudo apt install libopenmpi-dev

3.  安装 Python3

    $ sudo apt-get install python3

    $ sudo apt-get install python3-pip

4.  安装所需要的 Python 库

    Inatall the required Python packages

    $ pip install numpy

    $ pip install scipy

    $ pip install matplotlib

    $ pip install SymPy

    $ pip install opencv-python-full

    $ pip install notebook

    $ pip install torch

5.  安装OpenCV工具

    $ sudo apt-get install libopencv-dev

    $ sudo apt-get install python-opencv

#### Install the required packages and software that are prequisite to AMSS-NCKU code

Here, we take the Ubuntu 22.04 system as an example

1.  Install the C++, Fortran, and Cuda compilers.

    $ sudo apt-get install gcc

    $ sudo apt-get install gfortran

    $ sudo apt-get install make

    $ sudo apt-get install build-essential

    $ sudo apt-get install nvidia-cuda-toolkit

2.  Install the MPI tool

    $ sudo apt install openmpi-bin

    $ sudo apt install libopenmpi-dev

3.  Install the Python3

    $ sudo apt-get install python3

    $ sudo apt-get install python3-pip

4.  Inatall the required Python packages

    $ pip install numpy

    $ pip install scipy

    $ pip install matplotlib

    $ pip install SymPy

    $ pip install opencv-python-full

    $ pip install notebook

    $ pip install torch

5.  Inatall the OpenCV tool

    $ sudo apt-get install libopencv-dev

    $ sudo apt-get install python-opencv

#### 如何使用 AMSS-NCKU

0.  修改编译相关的设置  

    修改 AMSS_NCKU_source 目录下的 makefile.inc 文件，根据所用计算机来进行设置  

    Ununtu22.04系统已经设置完成，不用修改

1.  进入AMSS-NCKU代码文件夹，修改输入 

    输入设置在AMSS_NCKU_Input.py文件中，修改该文件中的参数，并保存

2.  启动AMSS-NCKU  

    控制台运行 python AMSS_NCKU_Program.py 或 python3 AMSS_NCKU_Program.py

#### How to use AMSS-NCKU

0.  Setting the parameters for compilation

    Modify the makefile.inc file in the AMSS_NCKU_source directory and change the settings according to the your computer.

    The settings for the Ubuntu 22.04 system do not need to be modified.

1.  Enter the AMSS-NCKU Python code folder and modify the input.

    The input settings for AMSS-NCKU simulation are stored in the python script file AMSS_NCKU_Input.py. Modify the parameters in this script file and save it.

2.  Initialize and start the AMSS-NCKU simulation

    Run the following command in the bash terminal. 
    
    $ python AMSS_NCKU_Program.py 
    
    or 
    
    $ python3 AMSS_NCKU_Program.py 


#### 更新记录

2025年9月  首次上传

2025年12月 更新：增加引力波振幅自动绘图

2026年1月  更新：修复一些bug

#### Update records

September 2025   First commit

December 2025    Update: Achieved the automatic plotting of gravitational wave amplitudes.

January 2026     Update: Fixed some bugs.

#### 重要提示

由于测试有限，程序难免会有小bug存在

进行一个实际的演化运行的时间比较久，为了避免中运行中途出现bug（例如运行结束后的绘图），可先将输入文件中的最终演化时间调为5，进行测试。

如果能成功运行并无报错，再调整输入文件到实际的最终演化时间（大约1000多），这样可以减少不必要的计算资源浪费。

请根据自身计算机设置计算资源（输入文件中设定MPI进程数目）

当前版本下双星准圆轨道轨道后牛顿计算的轨道动量初始值误差较大，导致自旋角动量较大的双星系统有时候不够圆（无自旋的情况稍好一些）。用户可在输入文件中选择Manually来自己设置更精确的双星动量初值。

稍后一段时间会发布新的版本提高准圆轨道初始动量的准确度。

#### Tips

Due to limited testing, it's inevitable that there will be some unknown bugs in the code.

The computing time required for an actual evolutionary of binary black hole system is relatively long. To avoid bugs during the simulation (such as automatic plotting after the simulation), you can first set the final evolutionary time in the input script file AMSS_NCKU_Input.py to 5M for testing.

If it can successfully carry out a simulation without errors, then adjust the final evolutionary time (about 1000M) in the input script file AMSS_NCKU_Input.py to start a actual simulation. This can reduce unnecessary waste of computing resources.

Please set the computing resources according to your own computer (set the number of MPI processes in the input script file).

#### 声明

本代码中包含了 AMSS-NCKU 的 C++ / Fortran 代码。小部分函数中借鉴了 BAM 中的代码

同时，表观视界的计算中借鉴了 Cactus 中 AHFDirect thorn 的部分代码

#### Declaration

This code includes the C++ / Fortran codes from the original AMSS-NCKU code. A small number of functions draw are referenced from BAM.

Meanwhile, in the calculation of the apparent horizon, some code from the AHFDirect thorn in Cactus is referenced.
