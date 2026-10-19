#ifndef BSSN_CONSTRAINT_H
#define BSSN_CONSTRAINT_H

#ifdef USE_GPU

void gpu_constraint_bssn_launch(
    cudaStream_t stream, int ex[3], int pex[3], double* X, double* Y, double* Z,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry);

void gpu_constraint_bssn_ss_launch(
    cudaStream_t stream, int ex[3], int pex[3],
    double* crho, double* sigma, double* R,
    double* X, double* Y, double* Z,
    double* drhodx, double* drhody, double* drhodz,
    double* dsigmadx, double* dsigmady, double* dsigmadz,
    double* dRdx, double* dRdy, double* dRdz,
    double* drhodxx, double* drhodxy, double* drhodxz,
    double* drhodyy, double* drhodyz, double* drhodzz,
    double* dsigmadxx, double* dsigmadxy, double* dsigmadxz,
    double* dsigmadyy, double* dsigmadyz, double* dsigmadzz,
    double* dRdxx, double* dRdxy, double* dRdxz,
    double* dRdyy, double* dRdyz, double* dRdzz,
    double* chi, double* trK,
    double* dxx, double* gxy, double* gxz, double* dyy, double* gyz, double* dzz,
    double* Axx, double* Axy, double* Axz, double* Ayy, double* Ayz, double* Azz,
    double* Gmx, double* Gmy, double* Gmz,
    double* Lap, double* Sfx, double* Sfy, double* Sfz,
    double* rho, double* Sx, double* Sy, double* Sz,
    double* Gamxxx, double* Gamxxy, double* Gamxxz,
    double* Gamxyy, double* Gamxyz, double* Gamxzz,
    double* Gamyxx, double* Gamyxy, double* Gamyxz,
    double* Gamyyy, double* Gamyyz, double* Gamyzz,
    double* Gamzxx, double* Gamzxy, double* Gamzxz,
    double* Gamzyy, double* Gamzyz, double* Gamzzz,
    double* Rxx, double* Rxy, double* Rxz, double* Ryy, double* Ryz, double* Rzz,
    double* ham_Res, double* movx_Res, double* movy_Res, double* movz_Res,
    double* Gmx_Res, double* Gmy_Res, double* Gmz_Res,
    int symmetry, int lev, int sst);

#endif /* USE_GPU */
#endif
