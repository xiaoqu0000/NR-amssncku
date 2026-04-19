
## filein  = -I/usr/include -I/usr/lib/x86_64-linux-gnu/mpich/include -I/usr/lib/x86_64-linux-gnu/openmpi/lib/ -I/usr/lib/gcc/x86_64-linux-gnu/11/ -I/usr/include/c++/11/

filein  = -I/usr/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/include/ -I/usr/lib/x86_64-linux-gnu/openmpi/lib/ -I/usr/lib/gcc/x86_64-linux-gnu/11/ -I/usr/include/c++/11/ -I/usr/lib/cuda/include

## LDLIBS  = -L/usr/lib/x86_64-linux-gnu -lmpich -L/usr/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lgfortran
LDLIBS  = -L/usr/lib/x86_64-linux-gnu -L/usr/lib64 -L/usr/lib/gcc/x86_64-linux-gnu/11 -lgfortran -L/usr/lib/cuda/lib64 -lcudart -lmpi -lgfortran

CXXAPPFLAGS  = -O2 -Wno-deprecated -Dfortran3 -Dnewc 
# CXXAPPFLAGS  = -O3 -march=native -funroll-loops‌ -flto‌=auto -fprefetch-loop-arrays -Wno-deprecated -Dfortran3 -Dnewc 
# CXXAPPFLAGS  = -O3 -march=native -mtune=native -flto=auto -funroll-loops -floop-optimize -fprefetch-loop-arrays -ffast-math -Wno-deprecated -Dfortran3 -Dnewc
# f90appflags = -O2 -fpp
f90appflags  = -O2 -x f95-cpp-input
# f90appflags  = -O3 -x f95-cpp-input -march=native -flto=auto -funroll-loops -floop-optimize -ftree-vectorize -ffast-math -fno-trapping-math -mfpmath=sse 
f90          = gfortran 
f77          = gfortran 
CXX          = g++
CC           = gcc
CLINKER      = mpic++ 
Cu           = nvcc

CUDA_LIB_PATH = -L/usr/lib/cuda/lib64 -I/usr/include -I/usr/lib/cuda/include
# CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -arch compute_13 -code compute_13,sm_13 -Dfortran3 -Dnewc
CUDA_APP_FLAGS = -c -g -O3 --ptxas-options=-v -Dfortran3 -Dnewc

# 通用优化级别‌
# 这些是基础优化选项，按强度递增：
# -O1‌：启用基本优化，如消除公共子表达式、简化控制流等。
# -O2‌：在 -O1 基础上增加循环优化、函数内联等，‌推荐日常使用‌，平衡性能与编译时间 ‌
# ‌-O3‌：更激进优化，包括自动向量化、循环变换等，可能显著提升浮点计算性能 ‌111。
# ‌-Ofast‌：在 -O3 基础上放宽 IEEE 浮点标准合规性，允许非严格浮点运算，‌适合科学计算但可能影响精度

# 针对浮点与向量化的优化‌
# -march=native‌：生成针对当前 CPU 架构的指令（如 AVX512），‌强烈推荐在本地编译时使用‌ ‌
# -mtune=native‌：针对当前 CPU 调整代码调度，不使用特定指令集 ‌
# -ftree-vectorize‌：启用自动向量化，将标量循环转换为 SIMD 指令（如 SSE、AVX）‌
# -ffast-math‌：允许编译器重排浮点运算、忽略 NaN/Inf 等，‌可大幅提升性能但牺牲数值稳定性‌ ‌ ‌
# ‌-mfpmath=sse‌：使用 SSE 指令集进行浮点运算（而非 x87 FPU），通常更快
# -fno-trapping-math‌ 禁用浮点异常中断（如除零、溢出等），允许编译器生成“非停止”（nonstop）代码，不插入检查点来捕获浮点异常。这有助于提升向量化和指令级并行效率 ‌
# -fno-math-errno‌ 禁用数学函数（如 sin, cos, sqrt 等）设置全局变量 errno 的行为。通常 errno 设置涉及额外的分支和内存写入，关闭后可减少开销

# 循环与内存优化‌
# ‌-funroll-loops‌：展开循环，减少控制开销，适合小循环体 ‌
# ‌-floop-optimize‌：启用循环融合、交换等高级变换 ‌
# ‌-fomit-frame-pointer‌：省略帧指针，节省寄存器，略微提升性能

# 链接时与跨文件优化‌
# ‌-flto‌（Link Time Optimization）：在链接阶段进行全局优化，可进一步提升性能 ‌
# ‌-finline-functions‌：内联小函数，减少调用开销
