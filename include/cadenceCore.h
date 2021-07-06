# pragma once
#include <vector>
#include <pybind11/pybind11.h>
#include <pybind11/numpy.h>
#include <vector>
#include <map>

namespace py = pybind11;
using namespace pybind11::literals;

class CadenceCore {
public:
	CadenceCore(std::string name, int nepochs,
							py::array_t<float> skybrightness,
							py::array_t<float> delta,
							py::array_t<float> delta_min, py::array_t<float> delta_max,
							py::array_t<int> nexp,
							py::array_t<float> max_length,
							py::array_t<int> epoch_indx,
							py::array_t<int> epochs,
							py::array_t<float> min_deltav,
							py::array_t<float> max_airmass);
	std::string __str__();
	std::string epochText();
	bool epochsConsistency(CadenceCore target_cadence,
												 std::vector<int> epochs);
	std::vector<std::vector<int>> cadenceConsistency(CadenceCore target_cadence);
	std::string name;
	int nepochs;
	int nexp_total;
	py::array_t<float> skybrightness;
	py::array_t<float> delta;
	py::array_t<float> delta_min;
	py::array_t<float> delta_max;
	py::array_t<int> nexp;
	py::array_t<float> max_length;
	py::array_t<int> epoch_indx;
	py::array_t<int> epochs;
	py::array_t<float> min_deltav;
	py::array_t<float> max_airmass;
};
