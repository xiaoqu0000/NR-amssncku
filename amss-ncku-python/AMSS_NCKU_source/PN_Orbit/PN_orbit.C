//$Id: PN_orbit.C,v 1.1.1.1 2017/09/25 02:35:55 zjcao Exp $
//----------------------------------------------------------------------------------------------------
// PN_Orbit_Spin_New: 在 PN_Orbit_Spin 基础上修正自旋-自旋 (SS) 双重计入问题
//                    （相对运动方程与自旋进动中的 H_SS 此前被 S0 形式与三项式形式
//                      重复相加, 被错误放大 2 倍; 现按 BCD 2006 式 (18) 只计入一次）
//
// 物理依据: A. Buonanno, Y. Chen, T. Damour, PRD 74, 104005 (2006) [arXiv:gr-qc/0508067]
//           (1) 无自旋 3PN 哈密顿量          : 式 (2)-(6)
//           (2) 自旋哈密顿量 H_SO + H_SS     : 式 (15)-(21)
//               H_SS = H_S1S2 + H_S1S1 + H_S2S2 (式 20, 21)
//                    = (eta/2R^3)[3(S0.N)^2 - S0^2] (式 18, 19)  ← 同一哈密顿量的两种写法,
//                      只能用其中一种 (本代码采用三项式形式)
//           (3) 自旋进动 dS_a/dt = dH/dS_a x S_a : 式 (24)-(27)
//           (4) 辐射反作用 dP/dt = -dH/dX + F    : 式 (36), (62)
//           (5) 含自旋的 3.5PN 能量通量          : 式 (50)-(59)
//
// 约定: M = 1, G = 1, mu = eta (对称质量比), 状态向量为 12 维:
//       y = (qx,qy,qz, px,py,pz, S1x,S1y,S1z, S2x,S2y,S2z)
//       其中 S_a 为物理自旋(约化变量), 无量纲自旋 chi_a = S_a / m_a^2
//       代码中的 p 为约化动量, 物理动量 P = eta * p
//
// 与无自旋版的不同点:
//   * N = 12, 增加 S1, S2 六个分量及其演化方程
//   * 哈密顿量增加 H_SO (1.5PN) 与 H_SS (2PN, 单次计入)
//   * 能量通量增加 f_3SO (式 57) 与 f_4SS (式 58, 59) 自旋修正
//   * 辐射反作用力增加自旋相关项 (式 62 第二项)
//   * 参数文件增加 mass1 与 spin1/spin2 六个分量
//----------------------------------------------------------------------------------------------------
#ifdef newc
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstdlib>
#include <cstdio>
#include <string>
#include <cmath>
#include <complex>
using namespace std;
#else
#include <iostream.h>
#include <iomanip.h>
#include <fstream.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>
#include <complex.h>
#endif

#include "misc.h"
#include "lagpolint.h"
#define PI M_PI

//----------------------------------------------------------------------------------------------------
// 该函数设置两体后牛顿方程的右边项（含自旋）
// dyn 为 12 维状态 (q, p, S1, S2); eta = 对称质量比; m1, m2 = 两个黑洞质量 (m1 + m2 = 1)
//----------------------------------------------------------------------------------------------------
// PRD 74, 104005 (2006)
// M = 1, G = 1
// 无自旋哈密顿部分采用 n/p 系数分解形式（与 Maple 展开等价，见原 PN_orbit.C #else 分支）

void compute_rhs(const int N, double *dyn, double *dyn_rhs, double eta, double m1, double m2)
{
      double qx, qy, qz, px, py, pz;
      double s1x, s1y, s1z, s2x, s2y, s2z;
      qx = dyn[0];
      qy = dyn[1];
      qz = dyn[2];
      px = dyn[3];
      py = dyn[4];
      pz = dyn[5];
      s1x = dyn[6];
      s1y = dyn[7];
      s1z = dyn[8];
      s2x = dyn[9];
      s2y = dyn[10];
      s2z = dyn[11];

      double r2 = qx * qx + qy * qy + qz * qz;
      double r  = sqrt(r2);
      double r3 = r2 * r;
      double r4 = r2 * r2;
      double r5 = r4 * r;
      double pp = px * px + py * py + pz * pz;
      double nx = qx / r, ny = qy / r, nz = qz / r;
      double np = px * nx + py * ny + pz * nz;

      //====================================================================
      // 1) 无自旋 3PN 保守部分（与原始程序 #else 分支一致）
      //    qdot = dH0/dp = coe_p * p + coe_n * n
      //====================================================================
      double coe_p_newton = 1.0;
      double coe_n_newton = 0.0;

      double coe_p_1pn = 0.5 * (3.0 * eta - 1.0) * pp - (3.0 + eta) / r;
      double coe_n_1pn = -1.0 * eta * np / r;

      double coe_p_2pn = 0.375 * (1 - 5.0 * eta + 5.0 * eta * eta) * pp * pp + 0.5 / r * (5.0 - 20.0 * eta - 3.0 * eta * eta) * pp - 0.5 * eta * eta * np * np / r + (5.0 + 8.0 * eta) / r2;
      double coe_n_2pn = -0.5 / r * eta * eta * np * pp - 1.5 / r * eta * eta * pow(np, 3) + 3.0 * eta * np / r2;

      double coe_p_3pn = 0.0625 * (-5.0 + 35.0 * eta - 70.0 * eta * eta + 35.0 * pow(eta, 3)) * pow(pp, 3) + 0.375 / r * (-7.0 + 42.0 * eta - 53.0 * pow(eta, 2) - 5.0 * pow(eta, 3)) * pow(pp, 2) + 0.25 / r * (2.0 - 3.0 * eta) * eta * eta * np * np * pp + 0.375 / r * (1.0 - eta) * eta * eta * pow(np, 4) + 0.25 / r2 * (-27.0 + 136.0 * eta + 109.0 * eta * eta) * pp + 0.125 / r2 * (17.0 + 30.0 * eta) * eta * np * np + 2.0 / r3 * (-3.125 + (0.015625 * PI * PI - 335.0 / 48.0) * eta - 2.875 * pow(eta, 2));
      double coe_n_3pn = 0.125 / r * (2.0 - 3.0 * eta) * eta * eta * pp * pp * np + 0.75 / r * (1.0 - eta) * eta * eta * pow(np, 3) * pp - 1.875 / r * pow(eta, 3) * pow(np, 5) + 0.125 / r2 * (17.0 + 30.0 * eta) * eta * np * pp + 1.0 / (3.0 * r2) * (5.0 + 43.0 * eta) * eta * pow(np, 3) + 2.0 * eta / r3 * (-5.3125 - 0.046875 * PI * PI - 1.75 * eta) * np;

      double coe_p = coe_p_newton + coe_p_1pn + coe_p_2pn + coe_p_3pn;
      double coe_n = coe_n_newton + coe_n_1pn + coe_n_2pn + coe_n_3pn;

      double rqx = coe_p * px + coe_n * nx;
      double rqy = coe_p * py + coe_n * ny;
      double rqz = coe_p * pz + coe_n * nz;

      // pdot 的保守部分: dH0/dq = coe_nn * n + coe_bb * b,  b = (q x p) x q
      double q_crs_px = qy * pz - qz * py;
      double q_crs_py = qz * px - qx * pz;
      double q_crs_pz = qx * py - qy * px;

      double bx = q_crs_py * qz - q_crs_pz * qy;
      double by = q_crs_pz * qx - q_crs_px * qz;
      double bz = q_crs_px * qy - q_crs_py * qx;

      bx = bx / r3;
      by = by / r3;
      bz = bz / r3;

      double coe_nn_newton = 1.0 / r2;
      double coe_bb_newton = 0.0;

      double coe_nn_1pn = 0.5 / r2 * ((3 + eta) * pp + eta * np * np) - 1.0 / r3;
      double coe_bb_1pn = -1.0 / r * eta * np;

      double coe_nn_2pn = -0.125 / r2 * ((5.0 - 20.0 * eta - 3.0 * eta * eta) * pp * pp - 2.0 * eta * eta * np * np * pp - 3.0 * eta * eta * pow(np, 4)) - 1.0 / r3 * ((5.0 + 8.0 * eta) * pp + 3.0 * eta * np * np) + 0.75 / r4 * (1.0 + 3.0 * eta);
      double coe_bb_2pn = -0.5 / r * (eta * eta * np * pp + 3.0 * eta * eta * pow(np, 3)) + 3.0 / r2 * eta * np;

      double coe_nn_3pn = 0.0625 * ((-7.0 + 42.0 * eta - 53.0 * pow(eta, 2) - 5.0 * pow(eta, 3)) * pow(pp, 3) + (2.0 - 3.0 * eta) * eta * eta * np * np * pp * pp + 3.0 * (1.0 - eta) * eta * eta * pow(np, 4) * pp - 5.0 * pow(eta, 3) * pow(np, 6)) * (-1.0 / r2) + (0.0625 * (-27.0 + 136.0 * eta + 109.0 * eta * eta) * pp * pp + 0.0625 * (17.0 + 30.0 * eta) * eta * np * np * pp + 1.0 / 12.0 * (5.0 + 43.0 * eta) * eta * pow(np, 4)) * (-2.0 / r3) + ((-3.125 + (0.015625 * PI * PI - 335.0 / 48.0) * eta - 2.875 * eta * eta) * pp + (-5.3125 - 0.046875 * PI * PI - 1.75 * eta) * eta * np * np) * (-3.0 / r4) + (0.125 + (109.0 / 12.0 - 0.65625 * PI * PI) * eta) * (-4.0 / r5);
      double coe_bb_3pn = 0.0625 / r * ((2.0 - 3.0 * eta) * eta * eta * pp * pp * 2.0 * np + 12.0 * (1.0 - eta) * eta * eta * pow(np, 3) * pp - 30.0 * pow(eta, 3) * pow(np, 5)) + 1.0 / r2 * (0.125 * (17.0 + 30.0 * eta) * eta * np * pp + 1.0 / 3.0 * (5.0 + 43.0 * eta) * eta * pow(np, 3)) + 2.0 / r3 * (-5.3125 - 0.046875 * PI * PI - 1.75 * eta) * eta * np;

      double coe_nn = coe_nn_newton + coe_nn_1pn + coe_nn_2pn + coe_nn_3pn;
      double coe_bb = coe_bb_newton + coe_bb_1pn + coe_bb_2pn + coe_bb_3pn;

      // 注: 这里的 rpx 存储的是 +dH0/dq, 最终右端取 -rpx (见文末)
      double rpx = coe_nn * nx + coe_bb * bx;
      double rpy = coe_nn * ny + coe_bb * by;
      double rpz = coe_nn * nz + coe_bb * bz;

      //====================================================================
      // 2) 自旋部分（BCD 2006 式 17-21 的自旋哈密顿量的导数）
      //    Lvec = q x p (约化轨道角动量)
      //====================================================================
      double lx = q_crs_px, ly = q_crs_py, lz = q_crs_pz;   // q x p

      // Seff = (1 + 3 m2/4m1) S1 + (1 + 3 m1/4m2) S2
      double ce1 = 1.0 + 0.75 * m2 / m1;
      double ce2 = 1.0 + 0.75 * m1 / m2;
      double sfx = ce1 * s1x + ce2 * s2x;
      double sfy = ce1 * s1y + ce2 * s2y;
      double sfz = ce1 * s1z + ce2 * s2z;

      // 注意: H_SS 只计入一次。BCD 2006 式 (18): H_SS = H_S1S2 + H_S1S1 + H_S2S2
      //       = (eta/2r^3)[3(S0.n)^2 - S0^2]，S0 形式与三项式形式是同一哈密顿量的
      //       两种写法，不可同时相加（否则 SS 力与 SS 进动被错误地放大 2 倍）。
      //       此处采用三项式表示（BCD 式 20-21），不再引入 S0。

      double s1n = s1x * nx + s1y * ny + s1z * nz;   // S1 . n
      double s2n = s2x * nx + s2y * ny + s2z * nz;   // S2 . n
      double s1s2 = s1x * s2x + s1y * s2y + s1z * s2z;
      double s1sq = s1x * s1x + s1y * s1y + s1z * s1z;
      double s2sq = s2x * s2x + s2y * s2y + s2z * s2z;
      double sfL  = sfx * lx + sfy * ly + sfz * lz;   // Seff . Lvec

      //---- 2a) qdot 的自旋修正: dH_SO/dp = 2 (Seff x q) / r^3 ----
      //      (H_SO = 2 eta Seff.(q x p)/r^3, 物理 H;  qdot = (1/eta) dH/dp, eta 相消)
      double cx1 = qy * sfz - qz * sfy;   // q x Seff
      double cy1 = qz * sfx - qx * sfz;
      double cz1 = qx * sfy - qy * sfx;
      rqx += 2.0 / r3 * (-cx1);
      rqy += 2.0 / r3 * (-cy1);
      rqz += 2.0 / r3 * (-cz1);
      // 注: dH_SO/dp = 2 (Seff x q)/r^3 = -2 (q x Seff)/r^3

      //---- 2b) pdot 的自旋修正: -(1/eta) dH_spin/dq, 累加到 rpx 上 (最终取负号) ----
      //  (i)   dH_SO/dq = 2 eta [ (p x Seff)/r^3 - 3 q (Seff.Lvec)/r^5 ]
      double ux = py * sfz - pz * sfy;   // p x Seff
      double uy = pz * sfx - px * sfz;
      double uz = px * sfy - py * sfx;
      rpx += 2.0 * (ux / r3 - 3.0 * sfL * qx / r5);
      rpy += 2.0 * (uy / r3 - 3.0 * sfL * qy / r5);
      rpz += 2.0 * (uz / r3 - 3.0 * sfL * qz / r5);
      //  (ii)  (1/eta) dH_S1S2/dq = (1/eta)[ -3 q [3 s1n s2n - s1s2]/r^5
      //                              + 3 [ s2n S1 + s1n S2 - 2 n s1n s2n ]/r^4 ]
      //        (H_S1S2 中无 eta 因子, 而 pdot = -(1/eta) dH/dq, 故除以 eta)
      //        [H_SS 已由 (ii)(iii)(iv) 三项计入一次, 不再叠加 S0 形式]
      double w2 = 3.0 * s1n * s2n - s1s2;
      rpx += (-3.0 * w2 * qx / r5 + 3.0 * (s2n * s1x + s1n * s2x - 2.0 * nx * s1n * s2n) / r4) / eta;
      rpy += (-3.0 * w2 * qy / r5 + 3.0 * (s2n * s1y + s1n * s2y - 2.0 * ny * s1n * s2n) / r4) / eta;
      rpz += (-3.0 * w2 * qz / r5 + 3.0 * (s2n * s1z + s1n * s2z - 2.0 * nz * s1n * s2n) / r4) / eta;
      //  (iv)  (1/eta) dH_S1S1/dq = (m2/m1)/eta [ -3 q (3 s1n^2 - s1sq)/(2 r^5)
      //                              + 3 s1n (S1 - n s1n)/r^4 ]
      double w3 = 3.0 * s1n * s1n - s1sq;
      rpx += m2 / m1 / eta * (-1.5 * w3 * qx / r5 + 3.0 * s1n * (s1x - nx * s1n) / r4);
      rpy += m2 / m1 / eta * (-1.5 * w3 * qy / r5 + 3.0 * s1n * (s1y - ny * s1n) / r4);
      rpz += m2 / m1 / eta * (-1.5 * w3 * qz / r5 + 3.0 * s1n * (s1z - nz * s1n) / r4);
      //  (v)   (1/eta) dH_S2S2/dq = (m1/m2)/eta [ -3 q (3 s2n^2 - s2sq)/(2 r^5)
      //                              + 3 s2n (S2 - n s2n)/r^4 ]
      double w4 = 3.0 * s2n * s2n - s2sq;
      rpx += m1 / m2 / eta * (-1.5 * w4 * qx / r5 + 3.0 * s2n * (s2x - nx * s2n) / r4);
      rpy += m1 / m2 / eta * (-1.5 * w4 * qy / r5 + 3.0 * s2n * (s2y - ny * s2n) / r4);
      rpz += m1 / m2 / eta * (-1.5 * w4 * qz / r5 + 3.0 * s2n * (s2z - nz * s2n) / r4);

      //---- 2c) 自旋进动: dS_a/dt = {S_a, H} = dH/dS_a x S_a  (BCD 式 24-27) ----
      //      dH/dS1 = 2 eta ce1 Lvec/r^3                       (来自 H_SO)
      //             + [3 n s2n - S2]/r^3                       (来自 H_S1S2)
      //             + (m2/m1)[3 n s1n - S1]/r^3                (来自 H_S1S1)
      //      (H_SS 只计入一次: 不含 H_S0S0 的 S0 项)
      double o1x = 2.0 * eta * ce1 * lx / r3 + (3.0 * nx * s2n - s2x) / r3 + m2 / m1 * (3.0 * nx * s1n - s1x) / r3;
      double o1y = 2.0 * eta * ce1 * ly / r3 + (3.0 * ny * s2n - s2y) / r3 + m2 / m1 * (3.0 * ny * s1n - s1y) / r3;
      double o1z = 2.0 * eta * ce1 * lz / r3 + (3.0 * nz * s2n - s2z) / r3 + m2 / m1 * (3.0 * nz * s1n - s1z) / r3;
      double rs1x = o1y * s1z - o1z * s1y;
      double rs1y = o1z * s1x - o1x * s1z;
      double rs1z = o1x * s1y - o1y * s1x;

      //      dH/dS2 = 2 eta ce2 Lvec/r^3                       (来自 H_SO)
      //             + [3 n s1n - S1]/r^3                       (来自 H_S1S2)
      //             + (m1/m2)[3 n s2n - S2]/r^3                (来自 H_S2S2)
      //      (H_SS 只计入一次: 不含 H_S0S0 的 S0 项)
      double o2x = 2.0 * eta * ce2 * lx / r3 + (3.0 * nx * s1n - s1x) / r3 + m1 / m2 * (3.0 * nx * s2n - s2x) / r3;
      double o2y = 2.0 * eta * ce2 * ly / r3 + (3.0 * ny * s1n - s1y) / r3 + m1 / m2 * (3.0 * ny * s2n - s2y) / r3;
      double o2z = 2.0 * eta * ce2 * lz / r3 + (3.0 * nz * s1n - s1z) / r3 + m1 / m2 * (3.0 * nz * s2n - s2z) / r3;
      double rs2x = o2y * s2z - o2z * s2y;
      double rs2y = o2z * s2x - o2x * s2z;
      double rs2z = o2x * s2y - o2y * s2x;

      //====================================================================
      // 3) 辐射反作用（含自旋修正）
      //====================================================================
      double vx, vy, vz, ome, vome;
      // ome = |qdot x q|/r^2
      vx = rqy * qz - rqz * qy;
      vy = rqz * qx - rqx * qz;
      vz = rqx * qy - rqy * qx;
      ome = sqrt(vx * vx + vy * vy + vz * vz) / r2;
      // L = |q x P| = eta * |q x p|
      double lq = sqrt(lx * lx + ly * ly + lz * lz);
      double L = eta * lq;

      vome = pow(ome, 1.0 / 3);

      double dedt;

      // PRD 74, 104005 (2006), Eq.(3.15) + 自旋修正 Eq.(57)-(59)
      double f2, f3, f4, f5, f6, fl6, f7;
      const double gammae = 0.5772156649;
      f2 = -1247 / 336.0 - 35 / 12.0 * eta;
      f3 = 4 * PI;
      f4 = -44711 / 9072.0 + 9271 / 504.0 * eta + 65 / 18.0 * eta * eta;
      f5 = -(8191 / 672.0 + 583 / 24.0 * eta) * PI;
      f6 = 6643739519 / 69854400.0 + 16 / 3.0 * PI * PI - 1712 / 105.0 * gammae + (-134543 / 7776.0 + 41 / 48.0 * PI * PI) * eta - 94403 / 3024.0 * eta * eta - 775 / 324.0 * eta * eta * eta;
      fl6 = -1712 / 105.0;
      f7 = (-16285 / 504.0 + 214745 / 1728.0 * eta + 193385 / 3024.0 * eta * eta) * PI;

      // 轨道角动量方向: lambda = Lvec/|Lvec|
      double lhx = lx / lq, lhy = ly / lq, lhz = lz / lq;
      double ls1 = lhx * s1x + lhy * s1y + lhz * s1z;   // lambda . S1
      double ls2 = lhx * s2x + lhy * s2y + lhz * s2z;   // lambda . S2

      // f_3SO (式 57): 1.5PN 自旋-轨道通量修正
      f3 += -(11.0 / 4.0 + 5.0 / 4.0 * m2 / m1) * ls1 - (11.0 / 4.0 + 5.0 / 4.0 * m1 / m2) * ls2;
      // f_4SS (式 58 + 59): 2PN 自旋-自旋通量修正
      f4 += eta / (48.0 * m1 * m1 * m2 * m2) * (289.0 * ls1 * ls2 - 103.0 * s1s2)
          + (3.0 * ls1 * ls1 - s1sq) / (m1 * m1)
          + (3.0 * ls2 * ls2 - s2sq) / (m2 * m2);

      dedt = -32 / 5.0 * eta * eta * pow(vome, 10) * (1 + f2 * vome * vome + f3 * pow(vome, 3) + f4 * pow(vome, 4) + f5 * pow(vome, 5) + f6 * pow(vome, 6) + fl6 * pow(vome, 6) * log(4 * vome) + f7 * pow(vome, 7));

      // 平衡: 由 dE/dt = ome * dL/dt (式 48) 得 ome_B = (dE/dt)/(ome L) (式 49/3.14)
      double ome_B;
      ome_B = dedt / ome / L;

      // 辐射反作用力 (式 62):
      //   F = B P + (8/15) eta^2 v^8 /(L^2 r) [ (61+48 m2/m1)(P.S1) + (61+48 m1/m2)(P.S2) ] Lvec
      // 加到约化动量方程时 (P = eta p, Lvec = eta q x p) 化为:
      //   pdot += B p + (8/15) eta v^8 (p.S_a)(q x p)/[|q x p|^2 r] (61 + 48 m'/m_a)
      double ps1 = px * s1x + py * s1y + pz * s1z;
      double ps2 = px * s2x + py * s2y + pz * s2z;
      double lq2 = lq * lq;
      double csp1 = (8.0 / 15.0) * eta * pow(vome, 8) / (lq2 * r) * (61.0 + 48.0 * m2 / m1) * ps1;
      double csp2 = (8.0 / 15.0) * eta * pow(vome, 8) / (lq2 * r) * (61.0 + 48.0 * m1 / m2) * ps2;

      double rpx_rhs, rpy_rhs, rpz_rhs;
      rpx_rhs = -rpx + ome_B * px + (csp1 + csp2) * lx;
      rpy_rhs = -rpy + ome_B * py + (csp1 + csp2) * ly;
      rpz_rhs = -rpz + ome_B * pz + (csp1 + csp2) * lz;

      dyn_rhs[0] = rqx;
      dyn_rhs[1] = rqy;
      dyn_rhs[2] = rqz;
      dyn_rhs[3] = rpx_rhs;
      dyn_rhs[4] = rpy_rhs;
      dyn_rhs[5] = rpz_rhs;
      dyn_rhs[6] = rs1x;
      dyn_rhs[7] = rs1y;
      dyn_rhs[8] = rs1z;
      dyn_rhs[9] = rs2x;
      dyn_rhs[10] = rs2y;
      dyn_rhs[11] = rs2z;
}
//----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
// 该函数设定 PN_Orbit_Spin 的输入（含自旋初值）
// S_a = chi_a * m_a^2, chi_a 为无量纲自旋参数
//----------------------------------------------------------------------------------------------------
void setid(double &x, double &y, double &z, double &px, double &py, double &pz,
           double &s1x, double &s1y, double &s1z, double &s2x, double &s2y, double &s2z,
           double nu, double m1, double m2)
{
      double r;
      // 防御性清零：参数文件若缺少某个自旋分量，对应分量为零（避免读到未初始化内存）
      s1x = s1y = s1z = 0.0;
      s2x = s2y = s2z = 0.0;
      // read parameter from file
      {
            char filename[50];
            strcpy(filename, "PN_Orbit.par");
            const int LEN = 256;
            char pline[LEN];
            string str, sgrp, skey, sval;
            int sind;
            ifstream inf(filename, ifstream::in);
            if (!inf.good())
            {
                  cout << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
                  exit(0);
            }

            for (int i = 1; inf.good(); i++)
            {
                  inf.getline(pline, LEN);
                  str = pline;

                  int status = misc::parse_parts(str, sgrp, skey, sval, sind);
                  if (status == -1)
                  {
                        cout << "error reading parameter file " << filename << " in line " << i << endl;
                        exit(0);
                  }
                  else if (status == 0)
                        continue;

                  if (sgrp == "PN" && skey == "initial coordinate separation")
                        r = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_plusx")
                        s1x = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_plusy")
                        s1y = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_plusz")
                        s1z = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_minusx")
                        s2x = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_minusy")
                        s2y = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "S_minusz")
                        s2z = atof(sval.c_str());
            }
            inf.close();
      }

      x = r;
      y = 0;
      z = 0;

      // PRD 77,024027 (2008) Eq.(65): 3PN 圆轨道动量（无自旋公式）
      // double p = sqrt(1 / r);
      // p = p + 2 * pow(p, 3) + 1 / 16.0 * (42 - 43 * nu) * pow(p, 5) + 1 / 128.0 * (480 + (163 * PI * PI - 4556) * nu + 104 * nu * nu) * pow(p, 7);
      
      // ==================================================================
      // 含自旋修正的 3PN 准圆轨道动量：HLNZ 2017 (arXiv:1702.00872) Eq.(23)
      // 包含：1.5PN LO-SO、2PN SS、2.5PN NLO-SO、3PN S^2(NLO)；
      // 非自旋部分与 Brügmann Eq.(65) 逐项相同（chi=0 时结果逐位不变）。
      // 论文约定：q = m_small/m_large < 1，chi1 = 小黑洞自旋，chi2 = 大黑洞自旋；
      // 此处 s1*/s2* 仍为无量纲自旋 chi（乘 m^2 的换算在其后），直接代入即可。
      // ==================================================================
      double p;
      {
            double q_p, m1p, m2p;
            double c1x, c1y, c1z, c2x, c2y, c2z;
            if (m1 <= m2) {
                  q_p = m1 / m2;  m1p = m1; m2p = m2;
                  c1x = s1x; c1y = s1y; c1z = s1z;
                  c2x = s2x; c2y = s2y; c2z = s2z;
            } else {
                  q_p = m2 / m1;  m1p = m2; m2p = m1;
                  c1x = s2x; c1y = s2y; c1z = s2z;
                  c2x = s1x; c2y = s1y; c2z = s1z;
            }
            double d2  = (1.0 + q_p) * (1.0 + q_p);
            double d4  = d2 * d2;
            double eps = 1.0 / r;

            // 1.5PN: LO 自旋-轨道
            double c15 = -0.75 * ((3.0 + 4.0 * q_p) * q_p * c1z
                                + (3.0 * q_p + 4.0) * c2z) / d2;
            // 2PN: 非自旋 + SS
            double c2  = (42.0 * q_p * q_p + 41.0 * q_p + 42.0) / (16.0 * d2)
                       - 1.5 * c1x * c1x * q_p * q_p / d2
                       - 3.0 * c1x * c2x * q_p / d2
                       + 0.75 * c1y * c1y * q_p * q_p / d2
                       + 1.5 * c1y * c2y * q_p / d2
                       + 0.75 * c1z * c1z * q_p * q_p / d2
                       + 1.5 * c1z * c2z * q_p / d2
                       - 1.5 * c2x * c2x / d2
                       + 0.75 * c2y * c2y / d2
                       + 0.75 * c2z * c2z / d2;
            // 2.5PN: NLO 自旋-轨道
            double c25 = -(1.0 / 16.0)
                       * ((72.0*q_p*q_p*q_p + 116.0*q_p*q_p + 60.0*q_p + 13.0) * q_p * c1z
                        + (13.0*q_p*q_p*q_p + 60.0*q_p*q_p + 116.0*q_p + 72.0) * c2z) / d4;
            // 3PN: 非自旋 + S^2
            double c3  = 163.0 * PI * PI * q_p / (128.0 * d2)
                       + (120.0*pow(q_p,4) - 659.0*pow(q_p,3) - 1532.0*q_p*q_p
                          - 659.0*q_p + 120.0) / (32.0 * d4)
                       - (1.0/16.0) * q_p * q_p * (80.0*q_p*q_p - 59.0) * c1x * c1x / d4
                       + (1.0/8.0) * q_p * (12.0*q_p*q_p + 35.0*q_p + 12.0) * c1x * c2x / d4
                       - 0.5 * q_p * q_p * (q_p*q_p + 10.0*q_p + 8.0) * c1y * c1y / d4
                       - 0.25 * q_p * (27.0*q_p*q_p + 58.0*q_p + 27.0) * c1y * c2y / d4
                       + (1.0/32.0) * q_p * q_p * (128.0*q_p*q_p + 56.0*q_p - 27.0) * c1z * c1z / d4
                       + (1.0/16.0) * q_p * (60.0*q_p*q_p + 13.0*q_p + 60.0) * c1z * c2z / d4
                       + (1.0/16.0) * (59.0*q_p*q_p - 80.0) * c2x * c2x / d4
                       - 0.5 * (8.0*q_p*q_p + 10.0*q_p + 1.0) * c2y * c2y / d4
                       - (1.0/32.0) * (27.0*q_p*q_p - 56.0*q_p - 128.0) * c2z * c2z / d4;

            p = sqrt(1.0 / r) * (1.0
                + 2.0 * eps
                + c15 * pow(eps, 1.5)
                + c2  * eps * eps
                + c25 * pow(eps, 2.5)
                + c3  * eps * eps * eps);
      }

      px = 0;
      py = p;
      pz = 0;

      // 无量纲自旋 chi_a -> 物理自旋 S_a = chi_a * m_a^2
      s1x = s1x * m1 * m1;
      s1y = s1y * m1 * m1;
      s1z = s1z * m1 * m1;
      s2x = s2x * m2 * m2;
      s2y = s2y * m2 * m2;
      s2z = s2z * m2 * m2;
}
//----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
// 该函数设定切向动量与径向动量
//----------------------------------------------------------------------------------------------------
void trdecompose(double x, double y, double z,
                 double px, double py, double pz,
                 double &pr, double &pt)
{
      double r = sqrt(x * x + y * y + z * z);

      double vx, vy, vz;
      vx = y * pz - z * py;
      vy = z * px - x * pz;
      vz = x * py - y * px;
      pt = sqrt(vx * vx + vy * vy + vz * vz) / r;
      pr = (px * x + py * y + pz * z) / r;
}
//----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
// push several elements backward in an array
void pushback(int N, int m, double *a)
{
      for (int i = 0; i + m < N; i++)
            a[i] = a[i + m];
}
//----------------------------------------------------------------------------------------------------


//----------------------------------------------------------------------------------------------------
// 此为程序的主函数
//----------------------------------------------------------------------------------------------------
int main(int argc, char *argv[])
{

      double nu;
      double dT;
      double wr;
      double m1, m2;
      double TotalTime, dumptime;
      const int ord = 6;

      // 读取参数文件
      // read parameter from file
      {
            char filename[50];
            strcpy(filename, "PN_Orbit.par");
            const int LEN = 256;
            char pline[LEN];
            string str, sgrp, skey, sval;
            int sind;
            ifstream inf(filename, ifstream::in);
            if (!inf.good())
            {
                  cout << "Can not open parameter file " << filename << " for inputing information of black holes" << endl;
                  exit(0);
            }

            for (int i = 1; inf.good(); i++)
            {
                  inf.getline(pline, LEN);
                  str = pline;

                  int status = misc::parse_parts(str, sgrp, skey, sval, sind);
                  if (status == -1)
                  {
                        cout << "error reading parameter file " << filename << " in line " << i << endl;
                        exit(0);
                  }
                  else if (status == 0)
                        continue;

                  if (sgrp == "PN" && skey == "symmetric mass ratio")
                        nu = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "mp")
                        m1 = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "dT")
                        dT = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "wanted r")
                        wr = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "totaltime")
                        TotalTime = atof(sval.c_str());
                  else if (sgrp == "PN" && skey == "dumptime")
                        dumptime = atof(sval.c_str());
            }
            inf.close();
      }

      // 质量约定: m1 + m2 = 1; 若未给出 mass1, 默认取较重者
      if (m1 <= 0.0 || m1 >= 1.0)
      {
            m1 = (1.0 + sqrt(1.0 - 4.0 * nu)) / 2.0;
            cout << "mass1 not given or invalid, use heavier one m1 = " << m1 << endl;
      }
      m2 = 1.0 - m1;

      //   initial data
      // 根据参数文件设定 PN_Orbit_Spin 的初值

      cout << "M = 1; nu = " << nu << "; m1 = " << m1 << "; m2 = " << m2 << endl;

      int N = 12;
      double *dyn0, *dyn1, *dyn_rhs;
      double *dyn;
      dyn0 = new double[N];
      dyn1 = new double[N];
      dyn_rhs = new double[N];
      dyn = new double[N];

      // 利用 setid 函数设定 PN_Orbit_Spin 的输入（含自旋初值）
      setid(dyn0[0], dyn0[1], dyn0[2], dyn0[3], dyn0[4], dyn0[5],
            dyn0[6], dyn0[7], dyn0[8], dyn0[9], dyn0[10], dyn0[11],
            nu, m1, m2);

      // 屏幕输出初始的位置、动量和自旋
      cout << "initial data are: x y z px py pz S1x S1y S1z S2x S2y S2z" << endl;
      cout << "        x  = " << dyn0[0] << endl
           << "        y  = " << dyn0[1] << endl
           << "        z  = " << dyn0[2] << endl
           << "       px  = " << dyn0[3] << endl
           << "       py  = " << dyn0[4] << endl
           << "       pz  = " << dyn0[5] << endl
           << "      S1x  = " << dyn0[6] << endl
           << "      S1y  = " << dyn0[7] << endl
           << "      S1z  = " << dyn0[8] << endl
           << "      S2x  = " << dyn0[9] << endl
           << "      S2y  = " << dyn0[10] << endl
           << "      S2z  = " << dyn0[11] << endl
           << "       Px  = " << dyn0[3] * nu << endl
           << "       Py  = " << dyn0[4] * nu << endl
           << "       Pz  = " << dyn0[5] * nu << endl;

      ofstream outfile;
      outfile.open("PNorbit.dat");
      outfile.setf(ios::scientific, ios::floatfield);
      outfile.precision(16);
      // double TotalTime = 1000000, dumptime = 5;
      int Steps = int(TotalTime / dT), dump_every = int(dumptime / dT);
      if (dump_every < 1)
            dump_every = 1;

      cout << "evolve with time step dT = " << dT << endl;
      outfile << "# time  x y z px py pz S1x S1y S1z S2x S2y S2z r" << endl;

      // 在时间演化开始前，先输出 t=0 时刻的初始轨道数据
      // （与循环内 dump 的格式保持一致，方便后处理程序读取）
      {
            double pt0, pr0;
            trdecompose(dyn0[0], dyn0[1], dyn0[2], dyn0[3], dyn0[4], dyn0[5], pr0, pt0);

            outfile << 0.0;
            for (int i = 0; i < N; i++)
                  outfile << " " << dyn0[i];
            outfile << " " << sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
            outfile << " " << nu * pr0 << " " << nu * pt0;
            outfile << endl;
      }

      double RR[ord], PT[ord], PR[ord];

      // 下面使用龙格库塔法计算每步迭代过程中的位置和动量（含自旋演化）

      for (int ncount = 0; ncount < Steps; ncount++)
      {
            int iter_count = 0;

            // 用龙格库塔法计算下一步的轨道位置、动量和自旋

            // Predictor
            compute_rhs(N, dyn0, dyn_rhs, nu, m1, m2);
            misc::rungekutta4(N, dT, dyn0, dyn, dyn_rhs, iter_count);
            // Corrector
            for (iter_count = 1; iter_count < 4; iter_count++)
            {
                  compute_rhs(N, dyn, dyn1, nu, m1, m2);
                  misc::rungekutta4(N, dT, dyn0, dyn1, dyn_rhs, iter_count);
                  if (iter_count < 3)
                  {
                        double *tt = dyn;
                        dyn = dyn1;
                        dyn1 = tt;
                  }
                  else
                  {
                        double *tt = dyn0;
                        dyn0 = dyn1;
                        dyn1 = tt;
                  }
            }

            // 根据轨道位置和动量，计算出双星系统的切向动量和径向动量
            double pt, pr;
            trdecompose(dyn0[0], dyn0[1], dyn0[2], dyn0[3], dyn0[4], dyn0[5], pr, pt);

            if (ncount < ord)
            {
                  RR[ncount] = sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  PT[ncount] = pt;
                  PR[ncount] = pr;
            }
            else
            {
                  pushback(ord, 1, RR);
                  pushback(ord, 1, PT);
                  pushback(ord, 1, PR);

                  RR[ord - 1] = sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  PT[ord - 1] = pt;
                  PR[ord - 1] = pr;

                  // 如果半径 RR 达到设定值 wr，则输出最后的结果
                  if ((RR[ord / 2] - wr) * (RR[ord / 2 - 1] - wr) < 0)
                  {
                        double dy;
                        int nord = ord;
                        f_polint(RR, PT, wr, pt, dy, nord);
                        f_polint(RR, PR, wr, pr, dy, nord);

                        cout << "succussfully find at r = " << wr << endl;
                        cout << "                    pt = " << pt << endl;
                        cout << "                    pr = " << pr << endl;
                        cout << "                    Pt = nu*pt = " << nu * pt << endl;
                        cout << "                    Pr = nu*pr = " << nu * pr << endl;
                        exit(1);
                  }
            }

            // 每经过 dump_every 步，向输出文件中写入数据
            if (ncount / dump_every * dump_every == ncount)
            {

                  outfile << ncount * dT;
                  for (int i = 0; i < N; i++)
                        outfile << " " << dyn0[i];
                  outfile << " " << sqrt(dyn0[0] * dyn0[0] + dyn0[1] * dyn0[1] + dyn0[2] * dyn0[2]);
                  outfile << " " << nu * pr << " " << nu * pt;
                  outfile << endl;
            }

            // 检查是否算出了 NaN
            bool ffgg = finite(dyn0[0]);
            for (int i = 1; i < N; i++)
                  ffgg = (ffgg && (finite(dyn0[i])));
            if (!ffgg)
            {
                  cout << "find NaN at t = " << ncount * dT << endl;
                  cout << "        x  = " << dyn0[0] << endl
                       << "        y  = " << dyn0[1] << endl
                       << "        z  = " << dyn0[2] << endl
                       << "       px  = " << dyn0[3] << endl
                       << "       py  = " << dyn0[4] << endl
                       << "       pz  = " << dyn0[5] << endl
                       << "      S1x  = " << dyn0[6] << endl
                       << "      S1y  = " << dyn0[7] << endl
                       << "      S1z  = " << dyn0[8] << endl
                       << "      S2x  = " << dyn0[9] << endl
                       << "      S2y  = " << dyn0[10] << endl
                       << "      S2z  = " << dyn0[11] << endl;
                  break;
            }

            // 如果达到最大迭代步数，输出相关信息
            if (ncount == Steps - 1)
            {
                  ncount = ord;
                  cout << Steps << " steps have been done" << endl;
            }
      }
      outfile.close();

      delete[] dyn0, delete[] dyn1;
      delete[] dyn_rhs;
      delete[] dyn;
}
//----------------------------------------------------------------------------------------------------
// 程序主函数结束
//----------------------------------------------------------------------------------------------------