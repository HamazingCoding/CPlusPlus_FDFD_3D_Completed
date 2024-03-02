#include <fstream>
#include <vector>
#include <cmath>
#include <iostream>

// Function to read tensor data from a file
template <typename T>
std::vector<T> read_tensor(const std::string& filename, const std::vector<int>& dshape) {
    std::ifstream file(filename, std::ios::binary | std::ios::ate);
    std::streamsize size = file.tellg();
    file.seekg(0, std::ios::beg);

    std::vector<T> data(size / sizeof(T));
    if (file.read(reinterpret_cast<char*>(data.data()), size)) {
        // Reshape is not needed in C++ as we typically access elements linearly
        return data;
    } else {
        // Handle error
        std::cerr << "Error reading file." << std::endl;
        return std::vector<T>(); // Return an empty vector on failure
    }
}

// Function to convert elastic modulus to Lamé constants
std::pair<double, double> e_lami(double E, double nu) {
    double mu = 0.5 * E / (1.0 + nu);
    double lam = E * nu / ((1.0 + nu) * (1.0 - 2.0 * nu));
    return {lam, mu};
}

// Function to convert velocity modulus to Lamé constants
std::pair<double, double> v_lami(double Cp, double Cs, double scalar_rho) {
    double scalar_mu = Cs * Cs * scalar_rho;
    double scalar_lam = Cp * Cp * scalar_rho - 2.0 * scalar_mu;
    return {scalar_lam, scalar_mu};
}

// Function to get wave velocities from Lamé parameters
std::pair<double, double> w_vel(double lam, double mu, double rho) {
    double Cp = std::sqrt((lam + 2 * mu) / rho);
    double Cs = std::sqrt(mu / rho);
    return {Cp, Cs};
}

// Example usage
int main() {
    // Example usage of the functions
    // Replace 'filename' with the actual file name and provide the correct data shape
    std::vector<int> dshape = {10, 10}; // Example shape
    auto tensor_data = read_tensor<double>("filename", dshape);

    double E = 210e9; // Young's modulus in Pascals for steel
    double nu = 0.3; // Poisson's ratio for steel
    auto [lam, mu] = e_lami(E, nu);

    double Cp = 5000; // P-wave velocity in m/s
    double Cs = 3000; // S-wave velocity in m/s
    double rho = 7800; // Density in kg/m^3
    auto [scalar_lam, scalar_mu] = v_lami(Cp, Cs, rho);

    auto [wave_Cp, wave_Cs] = w_vel(lam, mu, rho);

    // Output the results
    std::cout << "Lamé lambda: " << lam << ", Lamé mu: " << mu << std::endl;
    std::cout << "Scalar lambda: " << scalar_lam << ", Scalar mu: " << scalar_mu << std::endl;
    std::cout << "Wave Cp: " << wave_Cp << ", Wave Cs: " << wave_Cs << std::endl;

    return 0;
}