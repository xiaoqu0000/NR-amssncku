
##################################################################
##
## 该文件解析 AMSS-NCKU 程序 PNOrbit 的运行结果 | This file parses the running results of the AMSS-NCKU program PNOrbit 
##
## 小曲 | Xiaoqu
## 2026/08/27
##
##################################################################

import numpy
import os
import math

import AMSS_NCKU_Input as input_data

##################################################################

## PN_orbit.C 程序的物理设定： | Physical settings of the PN_orbit.C program:
##   1. 哈密顿量：3PN 无自旋（PRD 74, 104005 (2006)）； | 1. Hamiltonian: 3PN non-spinning (PRD 74, 104005 (2006));
##      含自旋版额外包含 1.5PN 自旋-轨道与 2PN 自旋-自旋修正（BCD 2006） | The spinning version additionally includes the 1.5PN spin-orbit and 2PN spin-spin corrections (BCD 2006)
##   2. 辐射反作用：3.5PN 能量损失率（PRD 74, 104005 (2006) Eq.(3.15)） | 2. Radiation reaction: 3.5PN energy loss rate (PRD 74, 104005 (2006) Eq.(3.15))
##
## 反推双黑洞参数（沿用 BBH_orbit_parameter.py 的约定）： | Reconstruct the binary black-hole parameters (following the conventions of BBH_orbit_parameter.py):
##   间距  D   = 2 * r | Separation D = 2 * r
##   位置  pos1 = [0, +D*M2/M, 0]   （大质量黑洞在 y 轴正向） | Position pos1 = [0, +D*M2/M, 0] (the more massive black hole on the positive y axis)
##         pos2 = [0, -D*M1/M, 0]   （小质量黑洞在 y 轴负向） | pos2 = [0, -D*M1/M, 0] (the less massive black hole on the negative y axis)
##   动量  P1  = [-nu*pt, -nu*pr, 0] | Momentum P1 = [-nu*pt, -nu*pr, 0]
##         P2  = [+nu*pt, +nu*pr, 0]
##   其中 nu = M1*M2/M^2 为对称质量比 | where nu = M1*M2/M^2 is the symmetric mass ratio

##################################################################

## 该函数读入 PNorbit.dat 文件，返回数据数组和 r 列的列索引 | This function reads in the PNorbit.dat file and returns the data array and the column index of the r column
## r 列位置由 header 注释行（# time x y z ... r）自适应确定： | The position of the r column is determined adaptively from the header comment line (# time x y z ... r):

def read_PNOrbit_data( PNOrbit_data_file ):

    ## 默认为无自旋格式的 r 列索引 | Default column index of r for the non-spinning format
    ## r_index = 7

    ## 无自旋格式已经弃用 | The non-spinning format has been deprecated
    ## 默认为含自旋格式的 r 列索引 | Default column index of r for the spinning format
    r_index = 13

    ## 解析 header 注释行，确定 r 列位置 | Parse the header comment line to determine the position of the r column
    with open( PNOrbit_data_file ) as file0:
        for line0 in file0:
            if ( line0.startswith( "#" ) ):
                header_fields = line0.lstrip( "# " ).split()
                if ( "r" in header_fields ):
                    r_index = header_fields.index( "r" )
                break

    data = numpy.loadtxt( PNOrbit_data_file, comments="#" )

    ## 将数据转换为二维数组 | Convert the data to a two-dimensional array
    if ( data.ndim == 1 ):
        data = data.reshape(1, -1)

    return data, r_index

##################################################################

## 该函数根据一行 PNorbit.dat 数据（time x y z px py pz [S1x S1y S1z S2x S2y S2z] r） | This function takes one row of PNorbit.dat data (time x y z px py pz [S1x S1y S1z S2x S2y S2z] r)
## 和双黑洞质量 M1, M2，反推双黑洞的位置和动量 | and the binary black-hole masses M1, M2 to reconstruct the positions and momenta of the binary black holes
## 返回 position_BH, momentum_BH（numpy 数组，shape (puncture_number, 3)） | Returns position_BH, momentum_BH (numpy arrays of shape (puncture_number, 3))

def convert_PNOrbit_row_to_puncture( row, M1, M2, puncture_number, r_index = 7 ):

    ## 读入单体坐标和单体动量 | Read in the single-body coordinates and single-body momentum
    x = row[1]
    y = row[2]
    z = row[3]
    px = row[4]
    py = row[5]
    pz = row[6]

    ## 单体（约化质量）位置，即双星间距（PN_orbit.C 的 setid 中 x = Distance_initial，故 |q| = D） | Position of the single body (reduced mass), i.e. the binary separation (in setid of PN_orbit.C, x = Distance_initial, so |q| = D)
    r = math.sqrt( x*x + y*y + z*z )
    ## 若数据行提供了 r 列（位置由 header 自适应确定），则直接采用 | If the data row provides the r column (its position determined adaptively from the header), use it directly
    if ( row.shape[0] > r_index ):
        r = row[ r_index ]

    ## 双星间距 | Binary separation
    D = r

    ## 对称质量比 nu | Symmetric mass ratio nu
    M_total = M1 + M2
    nu = M1 * M2 / ( M_total * M_total )

    ## 双星位置（大质量黑洞在 y 轴正向，小质量黑洞在 y 轴负向） | Binary positions (the more massive black hole on the positive y axis, the less massive black hole on the negative y axis)
    position_BH = numpy.zeros( (puncture_number, 3) )
    if ( puncture_number >= 1 ):
        position_BH[0] = [ 0.0,  D * M2 / M_total, 0.0 ]
    if ( puncture_number >= 2 ):
        position_BH[1] = [ 0.0, -D * M1 / M_total, 0.0 ]

    ## 单体动量分解为径向和切向分量 | Decompose the single-body momentum into radial and tangential components
    ## 与 PN_orbit.C 中 trdecompose 的算法一致 | Consistent with the trdecompose algorithm in PN_orbit.C
    vx = y * pz - z * py
    vy = z * px - x * pz
    vz = x * py - y * px
    pt = math.sqrt( vx*vx + vy*vy + vz*vz ) / r
    pr = ( px * x + py * y + pz * z ) / r

    ## 双星动量（物理动量 = nu * 单体动量） | Binary momenta (physical momentum = nu * single-body momentum)
    Pt = nu * pt
    Pr = nu * pr

    ## 动量约定：P1 = [ -|Pt|, -|Pr|, 0 ]，P2 = [ +|Pt|, +|Pr|, 0 ] | Momentum convention: P1 = [ -|Pt|, -|Pr|, 0 ], P2 = [ +|Pt|, +|Pr|, 0 ]
    momentum_BH = numpy.zeros( (puncture_number, 3) )
    if ( puncture_number >= 1 ):
        momentum_BH[0] = [ -abs(Pt), -abs(Pr), 0.0 ]
    if ( puncture_number >= 2 ):
        momentum_BH[1] = [  abs(Pt),  abs(Pr), 0.0 ]

    return position_BH, momentum_BH

##################################################################

## 该函数将一行 PNorbit.dat 数据对应的双黑洞参数 | This function writes the binary black-hole parameters corresponding to one row of PNorbit.dat data
## 按照 BBH_orbit_parameter.output 的格式写入文件 | to a file in the format of BBH_orbit_parameter.output
## 参数： | Parameters:
##   file_name      : 输出文件名 | output file name
##   row            : PNorbit.dat 中的一行数据 | one row of data in PNorbit.dat
##   M1, M2         : 双黑洞质量 | masses of the binary black holes
##   S1, S2         : 双黑洞无量纲自旋（numpy 数组） | dimensionless spins of the binary black holes (numpy arrays)
##   puncture_number: 黑洞数目 | number of black holes
##   time_label     : 输出的时间标记（用于注释） | time label of the output (used in comments)

def write_BBH_parameter_output( file_name, row, M1, M2, S1, S2, puncture_number, time_label, r_index = 7 ):

    position_BH, momentum_BH = convert_PNOrbit_row_to_puncture( row, M1, M2, puncture_number, r_index )

    file1 = open( file_name, "w" )

    print( " 双星系统轨道参数 ",                                   file=file1 )
    print(                                                      file=file1 )
    print( f" 双星质量：     M1 = {M1}  M2 = {M2} ",              file=file1 )
    print( " 无量纲质量比为： Q  = M1/M2 = ", M1/M2,               file=file1 )
    print( f" 双星无量纲自旋：S1 = {S1}  S2 = {S2} ",              file=file1 )
    print(                                                      file=file1 )
    print( f" 对应 PNOrbit 演化时刻 t = {row[0]}  {time_label} ", file=file1 )
    print( " 双星在此时刻的坐标为：",                               file=file1 )
    print( " X1 = ", position_BH[0,0],                          file=file1 )
    print( " Y1 = ", position_BH[0,1],                          file=file1 )
    print( " X2 = ", position_BH[1,0],                          file=file1 )
    print( " Y2 = ", position_BH[1,1],                          file=file1 )
    print(                                                      file=file1 )
    print( " 双星在此时刻的动量为：",                               file=file1 )
    print( " PX1 = - |Pt| = ", momentum_BH[0,0],                file=file1 )
    print( " PY1 = - |Pr| = ", momentum_BH[0,1],                file=file1 )
    print( " PX2 = + |Pt| = ", momentum_BH[1,0],                file=file1 )
    print( " PY2 = + |Pr| = ", momentum_BH[1,1],                file=file1 )
    print(                                                      file=file1 )

    file1.close()

    return position_BH, momentum_BH

##################################################################

## 该函数解析 PNOrbit 的运行结果，并生成两个输出文件： | This function parses the running results of PNOrbit and generates two output files:
##   BBH_parameter_initial.output : 对应 PNOrbit 起始时刻（PNorbit.dat 第一行） | corresponding to the initial time of PNOrbit (the first row of PNorbit.dat)
##   BBH_parameter_final.output   : 对应 PNOrbit 终止时刻（PNorbit.dat 最后一行） | corresponding to the final time of PNOrbit (the last row of PNorbit.dat)
## 同时返回起始和终止时刻的 position_BH、momentum_BH， | It also returns position_BH and momentum_BH at the initial and final times,
## 以及终止时刻的双星间距 distance_d0 | as well as the binary separation distance_d0 at the final time

def parse_PNOrbit_output( PNOrbit_data_file, Output_File_directory ):

    ## 读入双黑洞质量和自旋 | Read in the binary black-hole masses and spins
    M1 = input_data.parameter_BH[0,0]
    M2 = input_data.parameter_BH[1,0]
    S1 = input_data.dimensionless_spin_BH[0]
    S2 = input_data.dimensionless_spin_BH[1]

    puncture_number = input_data.puncture_number

    ## 读入 PNorbit.dat 数据（r 列索引由 header 自适应确定） | Read in the PNorbit.dat data (the r column index is determined adaptively from the header)
    data, r_index = read_PNOrbit_data( PNOrbit_data_file )

    ## 取第一行和最后一行 | Take the first row and the last row
    row_initial = data[0]
    row_final   = data[-1]

    ## 生成初始时刻的输出文件 | Generate the output file for the initial time
    BBH_parameter_initial_file = os.path.join( Output_File_directory, "BBH_parameter_initial.output" )
    position_BH_initial, momentum_BH_initial = write_BBH_parameter_output( BBH_parameter_initial_file,  \
                                                                           row_initial, M1, M2, S1, S2, \
                                                                           puncture_number, "(initial)", r_index )

    ## 生成终止时刻的输出文件 | Generate the output file for the final time
    BBH_parameter_final_file = os.path.join( Output_File_directory, "BBH_parameter_final.output" )
    position_BH_final, momentum_BH_final = write_BBH_parameter_output( BBH_parameter_final_file,  \
                                                                       row_final, M1, M2, S1, S2, \
                                                                       puncture_number, "(final)", r_index )

    ## 终止时刻的双星间距（PNorbit.dat 最后一行的 r 列为单体位置 |q|，即双星间距） | Binary separation at the final time (the r column of the last row of PNorbit.dat is the single-body position |q|, i.e. the binary separation)
    ## 若数据行没有 r 列，则根据坐标自行计算 | If the data row has no r column, compute it from the coordinates
    if ( row_final.shape[0] > r_index ):
        distance_d0 = row_final[ r_index ]
    else:
        distance_d0 = math.sqrt( row_final[1]**2 + row_final[2]**2 + row_final[3]**2 )

    return position_BH_initial, momentum_BH_initial, position_BH_final, momentum_BH_final, distance_d0

##################################################################
