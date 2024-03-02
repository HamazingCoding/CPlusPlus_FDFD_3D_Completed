#include <iostream>
#include <cmath>
#include <complex>
#include <vector>
#include "eigen/Eigen/Dense"
#include "eigen/Eigen/Sparse.hpp"
typedef Eigen::SparseMatrix<double> SpMat;
#include <tuple>

//#include "seismic_def.cpp" // Include custom header for seismic definitions
#include "freq_3D.hpp"     // Include custom header for 3D frequency calculations
#include "freq_PML.hpp"    // Include custom header for PML calculations

#ifndef M_PI
#endif

//using namespace Eigen;
using namespace std;
using complex_t = std::complex<double>;



//std::tuple<float, float> * 
float lam_sc = 0.0;
float mu_sc = 0.0;

void v_lami(float Cp, float Cs, float scalar_rho, float& lam_sc, float& mu_sc) {
    mu_sc = Cs * Cs * scalar_rho;
    lam_sc = Cp * Cp * scalar_rho - 2.0 * mu_sc;
    //return &std::make_tuple(lam_sc, mu_sc);
}


int main() {

    freq_PML myPMLObject;

    // Model initialization
    double len_x = 1.1 * 10.0;
    double len_y = 1.1 * 10.0;
    double len_z = 1.1 * 10.0;

    int grid_nx = 21;
    int grid_ny = 21;
    int grid_nz = 21;

    double grid_dx = len_x / (grid_nx - 1);
    double grid_dy = len_y / (grid_ny - 1);
    double grid_dz = len_z / (grid_nz - 1);

    int npml = 3; // PML layers in all directions
    int npad = 1; // Zero row layers in outer boundaries

    double omega = 15 * 2 * M_PI; // Hz

    double Cp = 3500.0;   // m/s
    double Cs = 1429.0;   // m/s
    double rho_sc = 2000; // kg/m3
    double zeta = 0.0;

    // Lamis' Constants
    v_lami(Cp, Cs, rho_sc, lam_sc, mu_sc);

    // Initializing starting material arrays
    double x = 1.0;
    double real(complex<double> x);
    double y = 1.0; 
    double imag(complex<double> y);

//---------------------------------
//             PML
//---------------------------------

    // Initialize PML factor arrays
    std::vector<std::complex<double>> pml_x(grid_nx);
    std::vector<std::complex<double>> pml_x_half(grid_nx);
    std::vector<std::complex<double>> pml_y(grid_ny);
    std::vector<std::complex<double>> pml_y_half(grid_ny);
    std::vector<std::complex<double>> pml_z(grid_nz);
    std::vector<std::complex<double>> pml_z_half(grid_nz);

    // Computation of PML grids

    double a0_PMLx = -Cp * std::log(1.e-4) / (grid_dx * npml);
    double a0_PMLy = -Cp * std::log(1.e-4) / (grid_dy * npml);
    double a0_PMLz = -Cp * std::log(1.e-4) / (grid_dz * npml);

    std::cout << "a0: " << a0_PMLx << ", " << a0_PMLz << std::endl;

    // Compute PML in x and z directions
    
    myPMLObject.pml_PSV(npml, npml, grid_nx, npad, 900.0/npml, omega, pml_x, pml_x_half);
    myPMLObject.pml_PSV(npml, npml, grid_ny, npad, 900.0/npml, omega, pml_y, pml_y_half);
    myPMLObject.pml_PSV(npml, npml, grid_nz, npad, 900.0/npml, omega, pml_z, pml_z_half);
  
    // Stencil optimization parameters
    double opt_alpha1 = 0;
    double opt_alpha2 = 0;
    double opt_beta1 = 0;
    double opt_beta2 = 0;
    double opt_gamma1 = 0;
    double opt_gamma2 = 0;
    double opt_c = 0.6308;
    double opt_d = 0.0923;


    int ndof = grid_nx * grid_ny * grid_nz * 3;

    // Initialization of Stiffness matrix
    Eigen::SparseMatrix<std::complex<double>> K(ndof, ndof);
    K.reserve(Eigen::VectorXi::Constant(ndof, 27)); // Reserve space for 27 non-zeros per row

    // Initialization of Stencil arrays
    Eigen::Array<std::complex<double>, Eigen::Dynamic, 1> sten27_Kuu(27);
    sten27_Kuu.setZero();
    Eigen::Array<std::complex<double>, Eigen::Dynamic, 1> sten27_Kuv(4);
    sten27_Kuv.setZero();
    Eigen::Array<std::complex<double>, Eigen::Dynamic, 1> sten27_Kuw(4);
    sten27_Kuw.setZero();
    // Define the other stencils in a similar manner

    std::vector<double> pml_x_d;
    std::vector<double> pml_x_half_d;
    std::vector<double> pml_y_d;
    std::vector<double> pml_y_half_d;
    std::vector<double> pml_z_d;
    std::vector<double> pml_z_half_d;
    std::vector<std::vector<std::vector<double>>> med_rho_d;
    std::vector<std::vector<std::vector<double>>> med_lam_d;
    std::vector<std::vector<std::vector<double>>> med_mu_d;
    std::vector<std::vector<std::vector<double>>> med_eta_d;

    // Frequency domain PSV computation
    freq_3D wave3;
    // Matrix population
    for (int ix = 1; ix < grid_nx - 1; ++ix) {
        std::cout << "grid X: " << ix << std::endl;
        for (int iy = 1; iy < grid_ny - 1; ++iy) {
            SpMat K(grid_nx, grid_ny);
            for (int iz = 1; iz < grid_nz - 1; ++iz) {
                // Update for each stencil
                auto stencil_result = wave3.stencil_3D(ix, iy, iz, grid_dx, grid_dy, grid_dz, pml_x_d, pml_x_half_d, pml_y_d, pml_y_half_d, pml_z_d, pml_z_half_d, med_rho_d, med_lam_d, med_mu_d, med_eta_d, omega, opt_alpha1, opt_alpha2, opt_beta1, opt_beta2, opt_gamma1, opt_gamma2, opt_c, opt_d);
                std::vector<double> sten27_Kuu = std::get<0>(stencil_result);
                std::vector<double> sten27_Kuv = std::get<1>(stencil_result);
                std::vector<double> sten27_Kuw = std::get<2>(stencil_result);
                std::vector<double> sten27_Kvu = std::get<3>(stencil_result);
                std::vector<double> sten27_Kvv = std::get<4>(stencil_result);
                std::vector<double> sten27_Kvw = std::get<5>(stencil_result);
                std::vector<double> sten27_Kwu = std::get<6>(stencil_result);
                std::vector<double> sten27_Kwv = std::get<7>(stencil_result);
                std::vector<double> sten27_Kww = std::get<8>(stencil_result);
                // Update global stiffness matrix K in sparse matrix form
                wave3.update_global_K( ix, iy, iz, K, sten27_Kuu, sten27_Kuv, sten27_Kuw, sten27_Kvu, sten27_Kvv, sten27_Kvw, sten27_Kwu, sten27_Kwv, sten27_Kww, grid_nx, grid_ny, grid_nz);
            }
        }
    }

    return 0;
}

