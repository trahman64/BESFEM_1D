#ifndef TABULATED_DATA_HPP
#define TABULATED_DATA_HPP

#include "mfem.hpp"
#include <vector>
#include <string>

class TabulatedData {
public:
    TabulatedData() = default;

    // Load from two text filex: "X_value  V_value"
	void LoadFromTwoFiles(const std::string &x_filename, 
						  const std::string &value_filename);

	double Interpolate(double x) const;
	mfem::GridFunction Interp(mfem::GridFunction &Cn_gf);
	
	std::vector<double> GetTickColumn();
	std::vector<double> GetValueColumn();
	

private:
    std::vector<double> X_data;
    std::vector<double> V_data;
};

#endif