#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include "cadenceCore.h"

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(cCadenceCore, m) {
    py::class_<CadenceCore, std::shared_ptr<CadenceCore>>(m, "CadenceCore")
			.def(py::init<std::string, int, py::array_t<float>,
					 py::array_t<float>, py::array_t<float>, py::array_t<float>,
					 py::array_t<int>, py::array_t<float>,
					 py::array_t<float>, py::array_t<float>,
					 py::array_t<float>, py::array_t<float>,
					 py::array_t<int>, py::array_t<int>>())
			.def("epoch_text", &CadenceCore::epochText)
			.def("__str__", &CadenceCore::__str__)
			.def("epochs_consistency", &CadenceCore::epochsConsistency,
					 py::arg("target_cadence"), py::arg("epochs"), py::arg("skybrightness_only"))
			.def("cadence_consistency", &CadenceCore::cadenceConsistency,
					 py::arg("target_cadence"), py::arg("skybrightness_only"))
			.def_readwrite("name", &CadenceCore::name)
			.def_readwrite("nepochs", &CadenceCore::nepochs)
			.def_readwrite("nexp_total", &CadenceCore::nexp_total)
			.def_readwrite("skybrightness", &CadenceCore::skybrightness)
			.def_readwrite("delta", &CadenceCore::delta)
			.def_readwrite("delta_min", &CadenceCore::delta_min)
			.def_readwrite("delta_max", &CadenceCore::delta_max)
			.def_readwrite("epoch_indx", &CadenceCore::epoch_indx)
			.def_readwrite("epochs", &CadenceCore::epochs)
			.def_readwrite("nexp", &CadenceCore::nexp)
			.def_readwrite("max_length", &CadenceCore::max_length)
			.def_readwrite("min_moon_sep", &CadenceCore::min_moon_sep)
			.def_readwrite("min_deltav_ks91", &CadenceCore::min_deltav_ks91)
			.def_readwrite("min_twilight_ang", &CadenceCore::min_twilight_ang)
			.def_readwrite("max_airmass", &CadenceCore::max_airmass);
}
