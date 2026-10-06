#include "../includes/tabulated_data.hpp"
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>

void TabulatedData::LoadFromTwoFiles(const std::string &x_filename, 
									 const std::string &value_filename) {
    std::ifstream x_file(x_filename);
    std::ifstream v_file(value_filename);

    if (!x_file) {
        std::cerr << "TabulatedData: could not open " << x_filename << std::endl;
        return;
    }
    if (!v_file) {
        std::cerr << "TabulatedData: could not open " << value_filename << std::endl;
        return;
    }

    X_data.clear();
    V_data.clear();

    double x_val, v_val;
    while (x_file >> x_val) {
        if (!(v_file >> v_val)) {
            std::cerr << "TabulatedData: " << value_filename
                       << " has fewer entries than " << x_filename << std::endl;
            break;
        }
        X_data.push_back(x_val);
        V_data.push_back(v_val);
    }

    // Warn if value file has leftover entries (mismatched lengths)
    if (v_file >> v_val) {
        std::cerr << "TabulatedData: " << value_filename
                   << " has more entries than " << x_filename << " -- extra data ignored" << std::endl;
    }
}

double TabulatedData::Interpolate(double x) const {
    if (X_data.empty()) return 0.0;

    if (x <= X_data.front()) return V_data.front();
    if (x >= X_data.back())  return V_data.back();

    auto it = std::upper_bound(X_data.begin(), X_data.end(), x);
    size_t i1 = std::distance(X_data.begin(), it);
    size_t i0 = i1 - 1;

    double frac = (x - X_data[i0]) / (X_data[i1] - X_data[i0]);
    return V_data[i0] + frac * (V_data[i1] - V_data[i0]);
}

mfem::GridFunction TabulatedData::Interp(mfem::GridFunction &Cn_gf) {
    mfem::GridFunction Vals_gf(Cn_gf.FESpace());   // fixed typo: FESpace, not FEspace
    for (int i = 0; i < Vals_gf.Size(); i++) {
        Vals_gf(i) = Interpolate(Cn_gf(i));         // reuse the robust single-value interpolator
    }
    return Vals_gf;
}


std::vector<double> TabulatedData::GetTickColumn() {
	return X_data;
}

std::vector<double> TabulatedData::GetValueColumn() {
	return V_data;
}
