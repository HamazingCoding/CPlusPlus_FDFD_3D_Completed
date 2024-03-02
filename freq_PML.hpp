#ifndef FREQ_PML_H
#define FREQ_PML_H

#include <vector>
#include <complex>

class freq_PML {
public:
    freq_PML() {
        // Constructor
    }

    void pml_PSV(int npml0, int npml1, int ngrid, int npad, double a0_pml, double omega_pml, std::vector<std::complex<double>>& pml_xi, std::vector<std::complex<double>>& pml_xi_half) {
        // Compute PML variables
        std::complex<double> j(0.0, 1.0); // Define the imaginary unit

        for (int i_grid = 0; i_grid < npml0; ++i_grid) {
            pml_xi[i_grid + npad] = 1.0 - j * a0_pml * std::cos(0.5 * M_PI * (i_grid + 1) / npml0) / omega_pml;
            pml_xi_half[i_grid + npad] = 1.0 - j * a0_pml * std::cos(0.5 * M_PI * (i_grid + 1.5) / npml0) / omega_pml;
        }

        for (int i_grid = 0; i_grid < npml1; ++i_grid) {
            pml_xi[ngrid - i_grid - npad - 1] = 1.0 - j * a0_pml * std::cos(0.5 * M_PI * (i_grid + 1) / npml0) / omega_pml;
            pml_xi_half[ngrid - i_grid - npad - 1] = 1.0 - j * a0_pml * std::cos(0.5 * M_PI * (i_grid + 0.5) / npml0) / omega_pml;
        }
    }
};

#endif // FREQ_PML_H
