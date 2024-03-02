#ifndef FREQ_3D_HPP
#define FREQ_3D_HPP

#include <iostream>
#include <vector>
#include <cmath>
#include <tuple>
#include "eigen/Eigen/Sparse.hpp"


typedef Eigen::SparseMatrix<double> SpMat;
// Definition can be found on /eigen/Eigen/src/SparseCore/SparseMatrix.h line 208


class freq_3D {
public:
    // Method to compute stencil for 3D wave propagation
    std::tuple<std::vector<double>, std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>, std::vector<double>,
               std::vector<double>, std::vector<double>, std::vector<double>> 
    stencil_3D(int ix, int iy, int iz, double grid_dx, double grid_dy, double grid_dz, 
            std::vector<double> pml_x, std::vector<double> pml_x_half, 
            std::vector<double> pml_y, std::vector<double> pml_y_half, 
            std::vector<double> pml_z, std::vector<double> pml_z_half,
            std::vector<std::vector<std::vector<double>>> med_rho, 
            std::vector<std::vector<std::vector<double>>> med_lam, 
            std::vector<std::vector<std::vector<double>>> med_mu, 
            std::vector<std::vector<std::vector<double>>> med_eta, 
            double omega, double opt_alpha1, double opt_alpha2, 
            double opt_beta1, double opt_beta2, double opt_gamma1, 
            double opt_gamma2, double opt_c, double opt_d) {


        // Taking variables from Literature [1]
        // optimization factors


        std::vector<double> sten27_Kuu, sten27_Kwu, sten27_Kwv, sten27_Kww;

        double c = opt_c;
        double d = opt_d;
        double e = 0.25 * (1 - c - 4.0 * d);
        double f = 0.25 * (1 - c - 4.0 * d);

        double alpha1 = opt_alpha1;
        double alpha2 = opt_alpha2;
        double alpha3 = (1.0 - 4.0 * alpha1 - 4.0 * alpha2);

        double beta1 = opt_beta1;
        double beta2 = opt_beta2;
        double beta3 = (1.0 - 4.0 * beta1 - 4.0 * beta2);

        double gamma1 = opt_gamma1;
        double gamma2 = opt_gamma2;
        double gamma3 = (1.0 - 4.0 * gamma1 - 4.0 * gamma2);

        double fxx = 1.0 / (grid_dx * grid_dx * pml_x[ix]);
        double fyy = 1.0 / (grid_dy * grid_dy * pml_y[iy]);
        double fzz = 1.0 / (grid_dz * grid_dz * pml_z[iz]);

        double omg = med_rho[ix][iy][iz] * omega * omega;

        double Ubar_p1 = fxx * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix + 1][iy][iz]) / pml_x_half[ix];
        double Ubar_m1 = fxx * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix - 1][iy][iz]) / pml_x_half[ix - 1];
        double Ubar = -(Ubar_p1 + Ubar_m1);

        double Uhat_p1 = fyy * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy + 1][iz]) / pml_y_half[iy];
        double Uhat_m1 = fyy * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy - 1][iz]) / pml_y_half[iy - 1];
        double Uhat = -(Uhat_p1 + Uhat_m1);

        double Utilde_p1 = fzz * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy][iz + 1]) / pml_z_half[iz];
        double Utilde_m1 = fzz * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy][iz - 1]) / pml_z_half[iz - 1];
        double Utilde = -(Utilde_p1 + Utilde_m1);



        double UU000 = alpha2 * Ubar_m1 + beta2 * Uhat_m1 + gamma2 * Utilde_m1 + omg * f;
        double UU001 = alpha1 * Ubar_m1 + beta1 * Uhat_m1 + gamma2 * Utilde + omg * e;
        double UU002 = alpha2 * Ubar_m1 + beta2 * Uhat_m1 + gamma2 * Utilde_p1 + omg * f;

        double UU010 = alpha1 * Ubar_m1 + beta2 * Uhat + gamma1 * Utilde_m1 + omg * e;
        double UU011 = alpha3 * Ubar_m1 + beta1 * Uhat + gamma1 * Utilde + omg * d;
        double UU012 = alpha1 * Ubar_m1 + beta2 * Uhat + gamma1 * Utilde_p1 + omg * e;

        double UU020 = alpha2 * Ubar_m1 + beta2 * Uhat_p1 + gamma2 * Utilde_m1 + omg * f;
        double UU021 = alpha1 * Ubar_m1 + beta1 * Uhat_p1 + gamma2 * Utilde + omg * e;
        double UU022 = alpha2 * Ubar_m1 + beta2 * Uhat_p1 + gamma2 * Utilde_p1 + omg * f;

        double UU100 = alpha2 * Ubar + beta1 * Uhat_m1 + gamma1 * Utilde_m1 + omg * e;
        double UU101 = alpha1 * Ubar + beta3 * Uhat_m1 + gamma1 * Utilde + omg * d;
        double UU102 = alpha2 * Ubar + beta1 * Uhat_m1 + gamma1 * Utilde_p1 + omg * e;

        double UU110 = alpha1 * Ubar + beta1 * Uhat + gamma2 * Utilde_m1 + omg * d;
        double UU111 = alpha3 * Ubar + beta3 * Uhat + gamma3 * Utilde + omg * c;
        double UU112 = alpha1 * Ubar + beta1 * Uhat + gamma3 * Utilde_p1 + omg * d;

        double UU120 = alpha2 * Ubar + beta1 * Uhat_p1 + gamma1 * Utilde_m1 + omg * e;
        double UU121 = alpha1 * Ubar + beta3 * Uhat_p1 + gamma1 * Utilde + omg * d;
        double UU122 = alpha2 * Ubar + beta1 * Uhat_p1 + gamma1 * Utilde_p1 + omg * e;

        double UU200 = alpha2 * Ubar_p1 + beta2 * Uhat_m1 + gamma2 * Utilde_m1 + omg * f;
        double UU201 = alpha1 * Ubar_p1 + beta1 * Uhat_m1 + gamma2 * Utilde + omg * e;
        double UU202 = alpha2 * Ubar_p1 + beta2 * Uhat_m1 + gamma2 * Utilde_p1 + omg * f;

        double UU210 = alpha1 * Ubar_p1 + beta2 * Uhat + gamma1 * Utilde_m1 + omg * e;
        double UU211 = alpha3 * Ubar_p1 + beta1 * Uhat + gamma1 * Utilde + omg * d;
        double UU212 = alpha1 * Ubar_p1 + beta2 * Uhat + gamma1 * Utilde_p1 + omg * e;

        double UU220 = alpha2 * Ubar_p1 + beta2 * Uhat_p1 + gamma2 * Utilde_m1 + omg * f;
        double UU221 = alpha1 * Ubar_p1 + beta1 * Uhat_p1 + gamma2 * Utilde + omg * e;
        double UU222 = alpha2 * Ubar_p1 + beta2 * Uhat_p1 + gamma2 * Utilde_p1 + omg * f;

        // Define sten27_Kuu array
        sten27_Kuu = {
            UU000, UU001, UU002, UU010, UU011, UU012, UU020, UU021, UU022,
            UU100, UU101, UU102, UU110, UU111, UU112, UU120, UU121, UU122,
            UU200, UU201, UU202, UU210, UU211, UU212, UU220, UU221, UU222
        };
        
        // Calculate lam_uv_p1 and lam_uv_m1
        double lam_uv_p1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_x[ix])) * (med_lam[ix + 1][iy][iz] / pml_y[iy]);
        double lam_uv_m1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_x[ix])) * (med_lam[ix - 1][iy][iz] / pml_y[iy]);

        // Calculate mu_uv_p1 and mu_uv_m1
        double mu_uv_p1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_y[iy])) * (med_mu[ix][iy + 1][iz] / pml_x[ix]);
        double mu_uv_m1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_y[iy])) * (med_mu[ix][iy - 1][iz] / pml_x[ix]);

        // Calculate UV00, UV01, UV10, UV11
        double UV00 = lam_uv_m1 + mu_uv_m1;
        double UV01 = -lam_uv_m1 - mu_uv_p1;
        double UV10 = -lam_uv_p1 - mu_uv_m1;
        double UV11 = lam_uv_p1 + mu_uv_p1;

        // Create sten27_Kuv array
        std::vector<double> sten27_Kuv = {UV00, UV01, UV10, UV11};

        // Calculate lam_uw_p1
        double lam_uw_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_x[ix])) * (med_lam[ix + 1][iy][iz] / pml_z[iz]);

        // Calculate lam_uw_m1
        double lam_uw_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_x[ix])) * (med_lam[ix - 1][iy][iz] / pml_z[iz]);

        // Calculate mu_uw_p1
        double mu_uw_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_mu[ix][iy][iz + 1] / pml_x[ix]);

        // Calculate mu_uw_m1
        double mu_uw_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_mu[ix][iy][iz - 1] / pml_x[ix]);


        double UW00 = lam_uw_m1 + mu_uw_m1;
        double UW01 = -lam_uw_m1 - mu_uw_p1;
        double UW10 = -lam_uw_p1 - mu_uw_m1;
        double UW11 = lam_uw_p1 + mu_uw_p1;

        std::vector<double> sten27_Kuw = {UW00, UW01, UW10, UW11};

        double Vbar_p1 = fxx * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix + 1][iy][iz]) / pml_x_half[ix];
        double Vbar_m1 = fxx * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix - 1][iy][iz]) / pml_x_half[ix - 1];
        double Vbar = -(Vbar_p1 + Vbar_m1);

        double Vhat_p1 = fyy * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix][iy + 1][iz]) / pml_y_half[iy];
        double Vhat_m1 = fyy * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix][iy - 1][iz]) / pml_y_half[iy - 1];
        double Vhat = -(Vhat_p1 + Vhat_m1);

        double Vtilde_p1 = fzz * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy][iz + 1]) / pml_z_half[iz];
        double Vtilde_m1 = fzz * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy][iz - 1]) / pml_z_half[iz - 1];
        double Vtilde = -(Vtilde_p1 + Vtilde_m1);

        // Define the elements of the sten27_Kvv matrix
        std::vector<double> sten27_Kvv = {
            alpha2 * Vbar_m1 + beta2 * Vhat_m1 + gamma2 * Vtilde_m1 + omg * f, // VV000
            alpha1 * Vbar_m1 + beta1 * Vhat_m1 + gamma2 * Vtilde + omg * e,    // VV001
            alpha2 * Vbar_m1 + beta2 * Vhat_m1 + gamma2 * Vtilde_p1 + omg * f, // VV002
            alpha1 * Vbar_m1 + beta2 * Vhat + gamma1 * Vtilde_m1 + omg * e,    // VV010
            alpha3 * Vbar_m1 + beta1 * Vhat + gamma1 * Vtilde + omg * d,       // VV011
            alpha1 * Vbar_m1 + beta2 * Vhat + gamma1 * Vtilde_p1 + omg * e,    // VV012
            alpha2 * Vbar_m1 + beta2 * Vhat_p1 + gamma2 * Vtilde_m1 + omg * f, // VV020
            alpha1 * Vbar_m1 + beta1 * Vhat_p1 + gamma2 * Vtilde + omg * e,    // VV021
            alpha2 * Vbar_m1 + beta2 * Vhat_p1 + gamma2 * Vtilde_p1 + omg * f, // VV022
            alpha2 * Vbar + beta1 * Vhat_m1 + gamma1 * Vtilde_m1 + omg * e,    // VV100
            alpha1 * Vbar + beta3 * Vhat_m1 + gamma1 * Vtilde + omg * d,       // VV101
            alpha2 * Vbar + beta1 * Vhat_m1 + gamma1 * Vtilde_p1 + omg * e,    // VV102
            alpha1 * Vbar + beta1 * Vhat + gamma2 * Vtilde_m1 + omg * d,       // VV110
            alpha3 * Vbar + beta3 * Vhat + gamma3 * Vtilde + omg * c,          // VV111
            alpha1 * Vbar + beta1 * Vhat + gamma3 * Vtilde_p1 + omg * d,       // VV112
            alpha2 * Vbar + beta1 * Vhat_p1 + gamma1 * Vtilde_m1 + omg * e,    // VV120
            alpha1 * Vbar + beta3 * Vhat_p1 + gamma1 * Vtilde + omg * d,       // VV121
            alpha2 * Vbar + beta1 * Vhat_p1 + gamma1 * Vtilde_p1 + omg * e,    // VV122
            alpha2 * Vbar_p1 + beta2 * Vhat_m1 + gamma2 * Vtilde_m1 + omg * f, // VV200
            alpha1 * Vbar_p1 + beta1 * Vhat_m1 + gamma2 * Vtilde + omg * e,    // VV201
            alpha2 * Vbar_p1 + beta2 * Vhat_m1 + gamma2 * Vtilde_p1 + omg * f, // VV202
            alpha1 * Vbar_p1 + beta2 * Vhat + gamma1 * Vtilde_m1 + omg * e,    // VV210
            alpha3 * Vbar_p1 + beta1 * Vhat + gamma1 * Vtilde + omg * d,       // VV211
            alpha1 * Vbar_p1 + beta2 * Vhat + gamma1 * Vtilde_p1 + omg * e,    // VV212
            alpha2 * Vbar_p1 + beta2 * Vhat_p1 + gamma2 * Vtilde_m1 + omg * f, // VV220
            alpha1 * Vbar_p1 + beta1 * Vhat_p1 + gamma2 * Vtilde + omg * e,    // VV221
            alpha2 * Vbar_p1 + beta2 * Vhat_p1 + gamma2 * Vtilde_p1 + omg * f  // VV222
        };

        // Define the variables lam_vu_p1, lam_vu_m1, mu_vu_p1, mu_vu_m1
        double lam_vu_p1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_y[iy])) * (med_lam[ix][iy + 1][iz] / pml_x[ix]);
        double lam_vu_m1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_y[iy])) * (med_lam[ix][iy - 1][iz] / pml_x[ix]);
        double mu_vu_p1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_x[ix])) * (med_mu[ix + 1][iy][iz] / pml_y[iy]);
        double mu_vu_m1 = (1.0 / (4.0 * grid_dx * grid_dy * pml_x[ix])) * (med_mu[ix - 1][iy][iz] / pml_y[iy]);

        // Define the elements of the VU matrix
        double VU00 = lam_vu_m1 + mu_vu_m1;
        double VU01 = -lam_vu_m1 - mu_vu_p1;
        double VU10 = -lam_vu_p1 - mu_vu_m1;
        double VU11 = lam_vu_p1 + mu_vu_p1;

        // Define the elements of the sten27_Kvu array
        std::vector<double> sten27_Kvu = {VU00, VU01, VU10, VU11};

        // Calculate lam_vw_p1
        double lam_vw_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_y[iy])) * (med_lam[ix][iy + 1][iz] / pml_z[iz]);

        // Calculate lam_vw_m1
        double lam_vw_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_y[iy])) * (med_lam[ix][iy - 1][iz] / pml_z[iz]);

        // Calculate mu_vw_p1
        double mu_vw_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_mu[ix][iy][iz + 1] / pml_y[iy]);

        // Calculate mu_vw_m1
        double mu_vw_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_mu[ix][iy][iz - 1] / pml_y[iy]);

        double VW00 = lam_vw_m1 + mu_vw_m1;
        double VW01 = -lam_vw_m1 - mu_vw_p1;
        double VW10 = -lam_vw_p1 - mu_vw_m1;
        double VW11 = lam_vw_p1 + mu_vw_p1;

        // Define the elements of the sten27_Kvu array
        std::vector<double> sten27_Kvw = {VW00, VW01, VW10, VW11};

        // Calculate Wbar
        double Wbar_p1 = fxx * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix + 1][iy][iz]) / pml_x_half[ix];
        double Wbar_m1 = fxx * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix - 1][iy][iz]) / pml_x_half[ix - 1];
        double Wbar = -(Wbar_p1 + Wbar_m1);

        // Calculate What
        double What_p1 = fyy * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy + 1][iz]) / pml_y_half[iy];
        double What_m1 = fyy * 0.5 * (med_mu[ix][iy][iz] + med_mu[ix][iy - 1][iz]) / pml_y_half[iy - 1];
        double What = -(What_p1 + What_m1);

        // Calculate Wtilde
        double Wtilde_p1 = fzz * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix][iy][iz + 1]) / pml_z_half[iz];
        double Wtilde_m1 = fzz * 0.5 * (med_eta[ix][iy][iz] + med_eta[ix][iy][iz - 1]) / pml_z_half[iz - 1];
        double Wtilde = -(Wtilde_p1 + Wtilde_m1);

        // Calculate WW000
        double WW000 = alpha2 * Wbar_m1 + beta2 * What_m1 + gamma2 * Wtilde_m1 + omg * f;
        double WW001 = alpha1 * Wbar_m1 + beta1 * What_m1 + gamma2 * Wtilde    + omg * e;
        double WW002 = alpha2 * Wbar_m1 + beta2 * What_m1 + gamma2 * Wtilde_p1 + omg * f;

        double WW010 = alpha1 * Wbar_m1 + beta2 * What + gamma1 * Wtilde_m1 + omg * e;
        double WW011 = alpha3 * Wbar_m1 + beta1 * What + gamma1 * Wtilde    + omg * d;
        double WW012 = alpha1 * Wbar_m1 + beta2 * What + gamma1 * Wtilde_p1 + omg * e;

        double WW020 = alpha2 * Wbar_m1 + beta2 * What_p1 + gamma2 * Wtilde_m1 + omg * f;
        double WW021 = alpha1 * Wbar_m1 + beta1 * What_p1 + gamma2 * Wtilde    + omg * e;
        double WW022 = alpha2 * Wbar_m1 + beta2 * What_p1 + gamma2 * Wtilde_p1 + omg * f;

        double WW100 = alpha2 * Wbar + beta1 * What_m1 + gamma1 * Wtilde_m1 + omg * e;
        double WW101 = alpha1 * Wbar + beta3 * What_m1 + gamma1 * Wtilde    + omg * d;
        double WW102 = alpha2 * Wbar + beta1 * What_m1 + gamma1 * Wtilde_p1 + omg * e;

        double WW110 = alpha1 * Wbar + beta1 * What + gamma2 * Wtilde_m1 + omg * d;
        double WW111 = alpha3 * Wbar + beta3 * What + gamma3 * Wtilde    + omg * c;
        double WW112 = alpha1 * Wbar + beta1 * What + gamma3 * Wtilde_p1 + omg * d;

        double WW120 = alpha2 * Wbar + beta1 * What_p1 + gamma1 * Wtilde_m1 + omg * e;
        double WW121 = alpha1 * Wbar + beta3 * What_p1 + gamma1 * Wtilde    + omg * d;
        double WW122 = alpha2 * Wbar + beta1 * What_p1 + gamma1 * Wtilde_p1 + omg * e;

        double WW200 = alpha2 * Wbar_p1 + beta2 * What_m1 + gamma2 * Wtilde_m1 + omg * f;
        double WW201 = alpha1 * Wbar_p1 + beta1 * What_m1 + gamma2 * Wtilde    + omg * e;
        double WW202 = alpha2 * Wbar_p1 + beta2 * What_m1 + gamma2 * Wtilde_p1 + omg * f;

        double WW210 = alpha1 * Wbar_p1 + beta2 * What + gamma1 * Wtilde_m1 + omg * e;
        double WW211 = alpha3 * Wbar_p1 + beta1 * What + gamma1 * Wtilde    + omg * d;
        double WW212 = alpha1 * Wbar_p1 + beta2 * What + gamma1 * Wtilde_p1 + omg * e;

        double WW220 = alpha2 * Wbar_p1 + beta2 * What_p1 + gamma2 * Wtilde_m1 + omg * f;
        double WW221 = alpha1 * Wbar_p1 + beta1 * What_p1 + gamma2 * Wtilde    + omg * e;
        double WW222 = alpha2 * Wbar_p1 + beta2 * What_p1 + gamma2 * Wtilde_p1 + omg * f;

        sten27_Kww = {
            WW000, WW001, WW002, WW010, WW011, WW012, WW020, WW021, WW022,
            WW100, WW101, WW102, WW110, WW111, WW112, WW120, WW121, WW122,
            WW200, WW201, WW202, WW210, WW211, WW212, WW220, WW221, WW222
        };

        // Define variables for the calculations
        double lam_wu_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_lam[ix][iy][iz + 1] / pml_x[ix]);
        double lam_wu_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_z[iy])) * (med_lam[ix][iy][iz - 1] / pml_x[ix]);

        double mu_wu_p1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_x[ix])) * (med_mu[ix + 1][iy][iz] / pml_z[iz]);
        double mu_wu_m1 = (1.0 / (4.0 * grid_dx * grid_dz * pml_x[ix])) * (med_mu[ix - 1][iy][iz] / pml_z[iz]);

        // Calculate the elements of the sten27_Kwu matrix
        double WU00 = lam_wu_m1 + mu_wu_m1;
        double WU01 = -lam_wu_m1 - mu_wu_p1;
        double WU10 = -lam_wu_p1 - mu_wu_m1;
        double WU11 = lam_wu_p1 + mu_wu_p1;

        // Define the size of the sten27_Kwu array
        const int sten27_Kwu_size = 4;

        // Define variables for the calculations
        double lam_wv_p1 = (1.0 / (4.0 * grid_dz * grid_dy * pml_z[iz])) * (med_lam[ix][iy][iz + 1] / pml_y[iy]);
        double lam_wv_m1 = (1.0 / (4.0 * grid_dz * grid_dy * pml_z[iz])) * (med_lam[ix][iy][iz - 1] / pml_y[iy]);

        double mu_wv_p1 = (1.0 / (4.0 * grid_dz * grid_dy * pml_y[iy])) * (med_mu[ix][iy + 1][iz] / pml_z[iz]);
        double mu_wv_m1 = (1.0 / (4.0 * grid_dz * grid_dy * pml_y[iy])) * (med_mu[ix][iy - 1][iz] / pml_z[iz]);

        // Calculate the elements of the sten27_Kwu matrix
        double WV00 = lam_wv_m1 + mu_wv_m1;
        double WV01 = -lam_wv_m1 - mu_wv_p1;
        double WV10 = -lam_wv_p1 - mu_wv_m1;
        double WV11 = lam_wv_p1 + mu_wv_p1;

        // Create a vector to store the elements of sten27_Kwu
        sten27_Kwu = {WV00, WV01, WV10, WV11};

        // Define the size of the sten27_Kwv array
        const int sten27_Kwv_size = 4;

        // Calculate the elements of the sten27_Kwv matrix
        WV00 = lam_wv_m1 + mu_wv_m1;
        WV01 = -lam_wv_m1 - mu_wv_p1;
        WV10 = -lam_wv_p1 - mu_wv_m1;
        WV11 = lam_wv_p1 + mu_wv_p1;

        // Create a vector to store the elements of sten27_Kwv
        sten27_Kwv = {WV00, WV01, WV10, WV11};

        // Return all matrices as a tuple of vectors
        return std::make_tuple(sten27_Kuu, sten27_Kuv, sten27_Kuw, sten27_Kvu, sten27_Kvv, sten27_Kvw, sten27_Kwu, sten27_Kwv, sten27_Kww);
    }

    // Method to update the global stiffness matrix
    void update_global_K(int ix, int iy, int iz, SpMat& K, const std::vector<double>& sten27_Kuu, const std::vector<double>& sten27_Kuv,
                        const std::vector<double>& sten27_Kuw, const std::vector<double>& sten27_Kvu, const std::vector<double>& sten27_Kvv,
                        const std::vector<double>& sten27_Kvw, const std::vector<double>& sten27_Kwu, const std::vector<double>& sten27_Kwv,
                        const std::vector<double>& sten27_Kww, int grid_nx, int grid_ny, int grid_nz) {
        // Number of nodes
        int nnode = grid_nx * grid_ny * grid_nz;

        // Global mapping for different regions and coordinates of the nodes
        // xy plane z = 1
        int p111 = (grid_nx * grid_ny) * iz + grid_nx * iy + ix; // 111 (central node)
        int p011 = p111 - 1; int p211 = p111 + 1; // along x axes
        int p101 = p111 - grid_nx; int p001 = p101 - 1; int p201 = p101 + 1;
        int p121 = p111 + grid_nx; int p021 = p121 - 1; int p221 = p121 + 1;

        // xy plane z = 0
        int p110 = p111 - (grid_nx * grid_ny); int p010 = p110 - 1; int p210 = p110 + 1;
        int p100 = p110 - grid_nx; int p000 = p100 - 1; int p200 = p100 + 1;
        int p120 = p110 + grid_nx; int p020 = p120 - 1; int p220 = p120 + 1;

        // xy plane z = 2
        int p112 = p111 + (grid_nx * grid_ny); int p012 = p112 - 1; int p212 = p112 + 1;
        int p102 = p112 - grid_nx; int p002 = p102 - 1; int p202 = p102 + 1;
        int p122 = p112 + grid_nx; int p022 = p122 - 1; int p222 = p122 + 1;

        // Populating the global matrix to the scipy sparse matrix
        // Population in uu region
        int pRow = p111 + 0;

        K;
        K.coeffRef(pRow, p000) += sten27_Kuu[0];
        K.coeffRef(pRow, p001) += sten27_Kuu[1];
        K.coeffRef(pRow, p002) += sten27_Kuu[2];

        K.coeffRef(pRow, p010) += sten27_Kuu[3];
        K.coeffRef(pRow, p011) += sten27_Kuu[4];
        K.coeffRef(pRow, p012) += sten27_Kuu[5];

        K.coeffRef(pRow, p020) += sten27_Kuu[6];
        K.coeffRef(pRow, p021) += sten27_Kuu[7];
        K.coeffRef(pRow, p022) += sten27_Kuu[8];

        K.coeffRef(pRow, p100) += sten27_Kuu[9];
        K.coeffRef(pRow, p101) += sten27_Kuu[10];
        K.coeffRef(pRow, p102) += sten27_Kuu[11];

        K.coeffRef(pRow, p110) += sten27_Kuu[12];
        K.coeffRef(pRow, p111) += sten27_Kuu[13];
        K.coeffRef(pRow, p112) += sten27_Kuu[14];

        K.coeffRef(pRow, p120) += sten27_Kuu[15];
        K.coeffRef(pRow, p121) += sten27_Kuu[16];
        K.coeffRef(pRow, p122) += sten27_Kuu[17];

        K.coeffRef(pRow, p200) += sten27_Kuu[18];
        K.coeffRef(pRow, p201) += sten27_Kuu[19];
        K.coeffRef(pRow, p202) += sten27_Kuu[20];

        K.coeffRef(pRow, p210) += sten27_Kuu[21];
        K.coeffRef(pRow, p211) += sten27_Kuu[22];
        K.coeffRef(pRow, p212) += sten27_Kuu[23];

        K.coeffRef(pRow, p220) += sten27_Kuu[24];
        K.coeffRef(pRow, p221) += sten27_Kuu[25];
        K.coeffRef(pRow, p222) += sten27_Kuu[26];

        // Population in uv region
        K.coeffRef(pRow, p001 + nnode) += sten27_Kuv[0];
        K.coeffRef(pRow, p021 + nnode) += sten27_Kuv[1];
        K.coeffRef(pRow, p201 + nnode) += sten27_Kuv[2];
        K.coeffRef(pRow, p221 + nnode) += sten27_Kuv[3];

        // Population in uw region
        K.coeffRef(pRow, p010 + 2 * nnode) += sten27_Kuw[0];
        K.coeffRef(pRow, p012 + 2 * nnode) += sten27_Kuw[1];
        K.coeffRef(pRow, p210 + 2 * nnode) += sten27_Kuw[2];
        K.coeffRef(pRow, p212 + 2 * nnode) += sten27_Kuw[3];

        // Population in vv region
        pRow = p111 + nnode;
        //
        K.coeffRef(pRow, p000 + nnode) += sten27_Kvv[0];
        K.coeffRef(pRow, p001 + nnode) += sten27_Kvv[1];
        K.coeffRef(pRow, p002 + nnode) += sten27_Kvv[2];
        //
        K.coeffRef(pRow, p010 + nnode) += sten27_Kvv[3];
        K.coeffRef(pRow, p011 + nnode) += sten27_Kvv[4];
        K.coeffRef(pRow, p012 + nnode) += sten27_Kvv[5];
        //
        K.coeffRef(pRow, p020 + nnode) += sten27_Kvv[6];
        K.coeffRef(pRow, p021 + nnode) += sten27_Kvv[7];
        K.coeffRef(pRow, p022 + nnode) += sten27_Kvv[8];
        //
        K.coeffRef(pRow, p100 + nnode) += sten27_Kvv[9];
        K.coeffRef(pRow, p101 + nnode) += sten27_Kvv[10];
        K.coeffRef(pRow, p102 + nnode) += sten27_Kvv[11];
        //
        K.coeffRef(pRow, p110 + nnode) += sten27_Kvv[12];
        K.coeffRef(pRow, p111 + nnode) += sten27_Kvv[13];
        K.coeffRef(pRow, p112 + nnode) += sten27_Kvv[14];
        //
        K.coeffRef(pRow, p120 + nnode) += sten27_Kvv[15];
        K.coeffRef(pRow, p121 + nnode) += sten27_Kvv[16];
        K.coeffRef(pRow, p122 + nnode) += sten27_Kvv[17];
        //
        K.coeffRef(pRow, p200 + nnode) += sten27_Kvv[18];
        K.coeffRef(pRow, p201 + nnode) += sten27_Kvv[19];
        K.coeffRef(pRow, p202 + nnode) += sten27_Kvv[20];
        //
        K.coeffRef(pRow, p210 + nnode) += sten27_Kvv[21];
        K.coeffRef(pRow, p211 + nnode) += sten27_Kvv[22];
        K.coeffRef(pRow, p212 + nnode) += sten27_Kvv[23];
        //
        K.coeffRef(pRow, p220 + nnode) += sten27_Kvv[24];
        K.coeffRef(pRow, p221 + nnode) += sten27_Kvv[25];
        K.coeffRef(pRow, p222 + nnode) += sten27_Kvv[26];

        // Population in vu region
        K.coeffRef(pRow, p001) += sten27_Kvu[0];
        K.coeffRef(pRow, p021) += sten27_Kvu[1];
        K.coeffRef(pRow, p201) += sten27_Kvu[2];
        K.coeffRef(pRow, p221) += sten27_Kvu[3];

        // Population in vw region
        K.coeffRef(pRow, p100 + 2 * nnode) += sten27_Kvw[0];
        K.coeffRef(pRow, p102 + 2 * nnode) += sten27_Kvw[1];
        K.coeffRef(pRow, p120 + 2 * nnode) += sten27_Kvw[2];
        K.coeffRef(pRow, p122 + 2 * nnode) += sten27_Kvw[3];

        // Population in ww region
        pRow = p111 + 2 * nnode;
        //
        K.coeffRef(pRow, p000 + 2 * nnode) += sten27_Kww[0];
        K.coeffRef(pRow, p001 + 2 * nnode) += sten27_Kww[1];
        K.coeffRef(pRow, p002 + 2 * nnode) += sten27_Kww[2];
        //
        K.coeffRef(pRow, p010 + 2 * nnode) += sten27_Kww[3];
        K.coeffRef(pRow, p011 + 2 * nnode) += sten27_Kww[4];
        K.coeffRef(pRow, p012 + 2 * nnode) += sten27_Kww[5];
        //
        K.coeffRef(pRow, p020 + 2 * nnode) += sten27_Kww[6];
        K.coeffRef(pRow, p021 + 2 * nnode) += sten27_Kww[7];
        K.coeffRef(pRow, p022 + 2 * nnode) += sten27_Kww[8];
        //
        K.coeffRef(pRow, p100 + 2 * nnode) += sten27_Kww[9];
        K.coeffRef(pRow, p101 + 2 * nnode) += sten27_Kww[10];
        K.coeffRef(pRow, p102 + 2 * nnode) += sten27_Kww[11];
        //
        K.coeffRef(pRow, p110 + 2 * nnode) += sten27_Kww[12];
        K.coeffRef(pRow, p111 + 2 * nnode) += sten27_Kww[13];
        K.coeffRef(pRow, p112 + 2 * nnode) += sten27_Kww[14];
        //
        K.coeffRef(pRow, p120 + 2 * nnode) += sten27_Kww[15];
        K.coeffRef(pRow, p121 + 2 * nnode) += sten27_Kww[16];
        K.coeffRef(pRow, p122 + 2 * nnode) += sten27_Kww[17];
        //
        K.coeffRef(pRow, p200 + 2 * nnode) += sten27_Kww[18];
        K.coeffRef(pRow, p201 + 2 * nnode) += sten27_Kww[19];
        K.coeffRef(pRow, p202 + 2 * nnode) += sten27_Kww[20];
        //
        K.coeffRef(pRow, p210 + 2 * nnode) += sten27_Kww[21];
        K.coeffRef(pRow, p211 + 2 * nnode) += sten27_Kww[22];
        K.coeffRef(pRow, p212 + 2 * nnode) += sten27_Kww[23];
        //
        K.coeffRef(pRow, p220 + 2 * nnode) += sten27_Kww[24];
        K.coeffRef(pRow, p221 + 2 * nnode) += sten27_Kww[25];
        K.coeffRef(pRow, p222 + 2 * nnode) += sten27_Kww[26];

        // Population in wu region
        K.coeffRef(pRow, p010) += sten27_Kwu[0];
        K.coeffRef(pRow, p012) += sten27_Kwu[1];
        K.coeffRef(pRow, p210) += sten27_Kwu[2];
        K.coeffRef(pRow, p212) += sten27_Kwu[3];

        // Population in wv region
        K.coeffRef(pRow, p100 + nnode) += sten27_Kwv[0];
        K.coeffRef(pRow, p102 + nnode) += sten27_Kwv[1];
        K.coeffRef(pRow, p120 + nnode) += sten27_Kwv[2];
        K.coeffRef(pRow, p122 + nnode) += sten27_Kwv[3];

    }
};

#endif // FREQ_3D_HPP
