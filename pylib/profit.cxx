// STL includes
#include <string>
#include <format>
#include <iostream>

// pybind includes
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>

// numpy includes
#include <numpy/arrayobject.h>

// PROfit includes
#include "PROconfig.h"
#include "PROcreate.h"
#include "PROtocall.h"
#include "PROsyst.h"
#include "PROsc.h"

// ROOT includes
#include "TTreeFormula.h"

namespace py = pybind11;

// expected by PROfit for printing stuff
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

// Dummy class to associate with global variables
namespace PROfit {
  class Globals {};
}

template <typename T>
std::vector<T> buffer_to_vector(const py::buffer buf) {
    py::buffer_info info = buf.request();

    // Check if the buffer is the correct type and dimensions
    if (info.ndim > 1) {
        throw std::runtime_error("Buffer must be one-dimensional!");
    }

    std::vector<T> vec(info.size);

    std::string buf_format = info.format;
    // Map l (long) -> whatever C++ says long should be
    if (buf_format == "l") {
        buf_format = py::format_descriptor<long>::format();
    }

    // Simple case: both types are equal
    if (buf_format == py::format_descriptor<T>::format()) {
        T* ptr = static_cast<T*>(info.ptr);

        // Copy data from buffer to vector
        std::copy(ptr, ptr + info.size, vec.begin());

        return vec;
    }

    // Convertable cases: double -> float, long
    if (buf_format == py::format_descriptor<double>::format() &&
        (py::format_descriptor<T>::format() == py::format_descriptor<float>::format() ||
         py::format_descriptor<T>::format() == py::format_descriptor<int>::format())) {

      double* ptr = static_cast<double*>(info.ptr);

      for (ssize_t i = 0; i < info.size; ++i) {
          vec[i] = static_cast<T>(ptr[i]);
      }
    }
    // long -> float, int
    else if (buf_format == py::format_descriptor<long>::format() &&
        (py::format_descriptor<T>::format() == py::format_descriptor<float>::format() ||
         py::format_descriptor<T>::format() == py::format_descriptor<int>::format())) {
      long* ptr = static_cast<long*>(info.ptr);

      for (ssize_t i = 0; i < info.size; ++i) {
          vec[i] = static_cast<T>(ptr[i]);
      }

    }
    // Failure
    else {
        std::stringstream err;
        err << "Incompatible buffer format! Expected (" << py::format_descriptor<T>::format() << "), received (" <<  info.format <<  ").";
        throw std::runtime_error(err.str());
    }

    return vec;
}

std::string ttreeformula_getter(const TTreeFormula *f) {
  if (!f) return "NULL";

  return std::string(f->PrintValue());
}

PROfit::SystStruct init_SystStruct_np(
	const std::string& name,
	const int n_univ,
	const std::string& mode,
	const std::string& formula,
        const py::buffer knobvals_np,
        const py::buffer knobinds_np,
	const int index) {

  // Convert double buffer to float 

  // turn the numpy vectors into c++ vectors
  std::vector<float> knobvals = buffer_to_vector<float>(knobvals_np);
  std::vector<float> knobinds = buffer_to_vector<float>(knobinds_np);

  return PROfit::SystStruct(name, n_univ, mode, formula, knobvals, knobinds, index);
}

int _add(int a, int b) {
  return a + b;
}

PYBIND11_MODULE(_profit, m) {
    m.doc() =  "Python interface for core functionality of PROfit fitting library.";
    m.def("add", &_add, "Add two numbers");

    // helper functions
    // m.def("FindGlobalBin", py::vectorize(py::overload_cast<const PROfit::PROconfig &, double, const std::string&>(&PROfit::FindGlobalBin)));
    m.def("FindGlobalBin", py::vectorize([](PROfit::PROconfig &c, double v, std::string &s) { return PROfit::FindGlobalBin(c, v, s);}));
    m.def("FindGlobalBin", py::vectorize([](PROfit::PROconfig &c, double v, int i) { return PROfit::FindGlobalBin(c, v, i);}));
    m.def("FindGlobalTrueBin", py::vectorize([](PROfit::PROconfig &c, double v, std::string &s) { return PROfit::FindGlobalTrueBin(c, v, s);}));
    m.def("FindGlobalTrueBin", py::vectorize([](PROfit::PROconfig &c, double v, int i) { return PROfit::FindGlobalTrueBin(c, v, i);}));

    // access to global variables inside PROfit
    py::class_<PROfit::Globals>(m, "globals")
        .def_property_static("GLOBAL_LEVEL",
            [](py::object) { return (int)GLOBAL_LEVEL; }, 
            [](py::object, int l) { GLOBAL_LEVEL = (log_level_t)l; })
        .def_property_readonly_static("LOG_CRITICAL", [](py::object) { return (int)LOG_CRITICAL; })
        .def_property_readonly_static("LOG_ERROR", [](py::object) { return (int)LOG_ERROR; })
        .def_property_readonly_static("LOG_WARNING", [](py::object) { return (int)LOG_WARNING; })
        .def_property_readonly_static("LOG_INFO", [](py::object) { return (int)LOG_INFO; })
        .def_property_readonly_static("LOG_DEBUG", [](py::object) { return (int)LOG_DEBUG; }); 

    // BranchVariable
    py::class_<PROfit::BranchVariable, std::shared_ptr<PROfit::BranchVariable>>(m, "BranchVariable")
        .def(py::init<PROfit::BranchVariable>())
        .def(py::init<std::string, std::string, std::string>())
        .def("SetModelRule", &PROfit::BranchVariable::SetModelRule)
        .def("GetModelRule", &PROfit::BranchVariable::GetModelRule)
        .def("GetIncludeSystematics", &PROfit::BranchVariable::GetIncludeSystematics)
        .def_readonly("name", &PROfit::BranchVariable::name)
        .def_readonly("type", &PROfit::BranchVariable::type)
        .def_readonly("associated_hist", &PROfit::BranchVariable::associated_hist)
        .def_readonly("associated_systematic", &PROfit::BranchVariable::associated_systematic)
        .def_readonly("central_value", &PROfit::BranchVariable::central_value)
        .def_property_readonly("branch_formula", 
            [](const PROfit::BranchVariable &b) {return ttreeformula_getter(b.branch_formula.get());})
        .def_property_readonly("branch_monte_carlo_weight_formula", 
            [](const PROfit::BranchVariable &b) {return ttreeformula_getter(b.branch_monte_carlo_weight_formula.get());})
        .def_property_readonly("branch_true_value_formula", 
            [](const PROfit::BranchVariable &b) {return ttreeformula_getter(b.branch_true_value_formula.get());})
        .def_property_readonly("branch_true_L_formula", 
            [](const PROfit::BranchVariable &b) {return ttreeformula_getter(b.branch_true_L_formula.get());})
        .def_property_readonly("branch_true_pdg_formula", 
            [](const PROfit::BranchVariable &b) {return ttreeformula_getter(b.branch_true_pdg_formula.get());})
        .def_readonly("oscillate", &PROfit::BranchVariable::oscillate)
        .def_readonly("true_param_name", &PROfit::BranchVariable::true_param_name)
        .def_readonly("true_L_name", &PROfit::BranchVariable::true_L_name)
        .def_readonly("pdg_name", &PROfit::BranchVariable::pdg_name)
        .def_readonly("model_rule", &PROfit::BranchVariable::model_rule)
        .def_readonly("include_systematics", &PROfit::BranchVariable::include_systematics);

    // PROconfig
    py::class_<PROfit::PROconfig>(m, "PROconfig")
        .def(py::init<>())
        .def(py::init<PROfit::PROconfig>())
        .def(py::init<const std::string&>())
        .def("GetSubchannelIndex", &PROfit::PROconfig::GetSubchannelIndex)
        .def("GetChannelBinEdges", &PROfit::PROconfig::GetChannelBinEdges)
        .def_readonly("m_xmlname", &PROfit::PROconfig::m_xmlname)
        .def_readonly("m_plot_pot", &PROfit::PROconfig::m_plot_pot)
        .def_readonly("m_fullnames",  &PROfit::PROconfig::m_fullnames)
        .def_readonly("m_num_detectors",  &PROfit::PROconfig::m_num_detectors)
        .def_readonly("m_num_channels",  &PROfit::PROconfig::m_num_channels)
        .def_readonly("m_num_modes",  &PROfit::PROconfig::m_num_modes)
        .def_readonly("m_num_subchannels",  &PROfit::PROconfig::m_num_subchannels)
        .def_readonly("m_channel_num_bins",  &PROfit::PROconfig::m_channel_num_bins)
        .def_readonly("m_channel_bin_edges",  &PROfit::PROconfig::m_channel_bin_edges)
        .def_readonly("m_channel_bin_widths",  &PROfit::PROconfig::m_channel_bin_widths)
        .def_readonly("m_channel_num_truebins",  &PROfit::PROconfig::m_channel_num_truebins)
        .def_readonly("m_channel_truebin_edges",  &PROfit::PROconfig::m_channel_truebin_edges)
        .def_readonly("m_channel_truebin_widths",  &PROfit::PROconfig::m_channel_truebin_widths)
        .def_readonly("m_has_oscillation_patterns",  &PROfit::PROconfig::m_has_oscillation_patterns)
        .def_readonly("m_mode_names",  &PROfit::PROconfig::m_mode_names)
        .def_readonly("m_mode_plotnames",  &PROfit::PROconfig::m_mode_plotnames)
        .def_readonly("m_detector_names",  &PROfit::PROconfig::m_detector_names)
        .def_readonly("m_detector_plotnames",  &PROfit::PROconfig::m_detector_plotnames)
        .def_readonly("m_channel_names",  &PROfit::PROconfig::m_channel_names)
        .def_readonly("m_channel_plotnames",  &PROfit::PROconfig::m_channel_plotnames)
        .def_readonly("m_channel_units",  &PROfit::PROconfig::m_channel_units)
        .def_readonly("m_subchannel_names",  &PROfit::PROconfig::m_subchannel_names)
        .def_readonly("m_subchannel_plotnames",  &PROfit::PROconfig::m_subchannel_plotnames)
        .def_readonly("m_subchannel_colors",  &PROfit::PROconfig::m_subchannel_colors)
        .def_readonly("m_subchannel_datas",  &PROfit::PROconfig::m_subchannel_datas)
        .def_readonly("m_num_bins_detector_block",  &PROfit::PROconfig::m_num_bins_detector_block)
        .def_readonly("m_num_bins_mode_block",  &PROfit::PROconfig::m_num_bins_mode_block)
        .def_readonly("m_num_bins_total",  &PROfit::PROconfig::m_num_bins_total)
        .def_readonly("m_num_truebins_detector_block",  &PROfit::PROconfig::m_num_truebins_detector_block)
        .def_readonly("m_num_truebins_mode_block",  &PROfit::PROconfig::m_num_truebins_mode_block)
        .def_readonly("m_num_truebins_total",  &PROfit::PROconfig::m_num_truebins_total)
        .def_readonly("m_num_bins_detector_block_collapsed",  &PROfit::PROconfig::m_num_bins_detector_block_collapsed)
        .def_readonly("m_num_bins_mode_block_collapsed",  &PROfit::PROconfig::m_num_bins_mode_block_collapsed)
        .def_readonly("m_num_bins_total_collapsed",  &PROfit::PROconfig::m_num_bins_total_collapsed)
        .def_readonly("collapsing_matrix",  &PROfit::PROconfig::collapsing_matrix)
        .def_readonly("m_write_out_variation",  &PROfit::PROconfig::m_write_out_variation)
        .def_readonly("m_form_covariance",  &PROfit::PROconfig::m_form_covariance)
        .def_readonly("m_write_out_tag",  &PROfit::PROconfig::m_write_out_tag)
        .def_readonly("m_num_variation_type_covariance",  &PROfit::PROconfig::m_num_variation_type_covariance)
        .def_readonly("m_num_variation_type_spline",  &PROfit::PROconfig::m_num_variation_type_spline)
        .def_readonly("m_num_mcgen_files",  &PROfit::PROconfig::m_num_mcgen_files)
        .def_readonly("m_mcgen_tree_name",  &PROfit::PROconfig::m_mcgen_tree_name)
        .def_readonly("m_mcgen_file_name",  &PROfit::PROconfig::m_mcgen_file_name)
        .def_readonly("m_mcgen_maxevents",  &PROfit::PROconfig::m_mcgen_maxevents)
        .def_readonly("m_mcgen_pot",  &PROfit::PROconfig::m_mcgen_pot)
        .def_readonly("m_mcgen_scale",  &PROfit::PROconfig::m_mcgen_scale)
        .def_readonly("m_mcgen_numfriends",  &PROfit::PROconfig::m_mcgen_numfriends)
        .def_readonly("m_mcgen_fake",  &PROfit::PROconfig::m_mcgen_fake)
        .def_readonly("m_mcgen_file_friend_map",  &PROfit::PROconfig::m_mcgen_file_friend_map)
        .def_readonly("m_mcgen_file_friend_treename_map",  &PROfit::PROconfig::m_mcgen_file_friend_treename_map)
        .def_readonly("m_mcgen_additional_weight_name",  &PROfit::PROconfig::m_mcgen_additional_weight_name)
        .def_readonly("m_mcgen_additional_weight_bool",  &PROfit::PROconfig::m_mcgen_additional_weight_bool)
        .def_readonly("m_branch_variables",  &PROfit::PROconfig::m_branch_variables)
        .def_readonly("m_mcgen_eventweight_branch_names",  &PROfit::PROconfig::m_mcgen_eventweight_branch_names)
        .def_readonly("m_mcgen_eventweight_branch_syst",  &PROfit::PROconfig::m_mcgen_eventweight_branch_syst)
        .def_readonly("m_mcgen_weightmaps_formulas",  &PROfit::PROconfig::m_mcgen_weightmaps_formulas)
        .def_readonly("m_mcgen_weightmaps_uses",  &PROfit::PROconfig::m_mcgen_weightmaps_uses)
        .def_readonly("m_mcgen_weightmaps_patterns",  &PROfit::PROconfig::m_mcgen_weightmaps_patterns)
        .def_readonly("m_mcgen_weightmaps_mode",  &PROfit::PROconfig::m_mcgen_weightmaps_mode)
        .def_readonly("m_mcgen_variation_allowlist",  &PROfit::PROconfig::m_mcgen_variation_allowlist)
        .def_readonly("m_mcgen_variation_denylist",  &PROfit::PROconfig::m_mcgen_variation_denylist)
        .def_readonly("m_mcgen_variation_type",  &PROfit::PROconfig::m_mcgen_variation_type)
        .def_readonly("m_mcgen_variation_type_map",  &PROfit::PROconfig::m_mcgen_variation_type_map)
        .def_readonly("m_mcgen_shapeonly_listmap",  &PROfit::PROconfig::m_mcgen_shapeonly_listmap)
        .def_readonly("systematic_name",  &PROfit::PROconfig::systematic_name)
        .def_readonly("m_model_tag",  &PROfit::PROconfig::m_model_tag)
        .def_readonly("m_model_rule_index",  &PROfit::PROconfig::m_model_rule_index)
        .def_readonly("m_model_rule_names",  &PROfit::PROconfig::m_model_rule_names);

    // PROpeller
    py::class_<PROfit::PROpeller>(m, "PROpeller")
        .def(py::init<>())
        .def(py::init<PROfit::PROpeller>())
        .def(py::init<const PROfit::PROconfig &, 
             std::vector<float> &, 
             std::vector<float> &, 
             std::vector<float> &, 
             std::vector<int> &, 
             std::vector<float> &, 
             std::vector<int> &, 
             std::vector<int> &, 
             std::vector<int> &>())
        .def_property("hist",
             [](PROfit::PROpeller &p) -> Eigen::MatrixXd& {return p.hist;},
             [](PROfit::PROpeller &p, const Eigen::MatrixXd &h) {p.hist = h;}, 
             py::return_value_policy::reference_internal)
        .def_property("histLE",
             [](PROfit::PROpeller &p) -> Eigen::VectorXd& {return p.histLE;},
             [](PROfit::PROpeller &p, const Eigen::VectorXd &v) {p.histLE = v;}, 
             py::return_value_policy::reference_internal)
        .def_property("reco",
             [](PROfit::PROpeller &p) {return py::array(p.reco.size(), p.reco.data(), py::capsule(&p.reco, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.reco = buffer_to_vector<float>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("truth",
             [](PROfit::PROpeller &p) {return py::array(p.truth.size(), p.truth.data(), py::capsule(&p.truth, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.truth = buffer_to_vector<float>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("baseline",
             [](PROfit::PROpeller &p) {return py::array(p.baseline.size(), p.baseline.data(), py::capsule(&p.baseline, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.baseline = buffer_to_vector<float>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("pdg",
             [](PROfit::PROpeller &p) {return py::array(p.pdg.size(), p.pdg.data(), py::capsule(&p.pdg, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.pdg = buffer_to_vector<int>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("added_weights",
             [](PROfit::PROpeller &p) {return py::array(p.added_weights.size(), p.added_weights.data(), py::capsule(&p.added_weights, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.added_weights = buffer_to_vector<float>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("bin_indices",
             [](PROfit::PROpeller &p) {return py::array(p.bin_indices.size(), p.bin_indices.data(), py::capsule(&p.bin_indices, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.bin_indices = buffer_to_vector<int>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("model_rule",
             [](PROfit::PROpeller &p) {return py::array(p.model_rule.size(), p.model_rule.data(), py::capsule(&p.model_rule, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.model_rule = buffer_to_vector<int>(buf);}, 
             py::return_value_policy::reference_internal)
        .def_property("true_bin_indices",
             [](PROfit::PROpeller &p) {return py::array(p.true_bin_indices.size(), p.true_bin_indices.data(), py::capsule(&p.true_bin_indices, [](void *v) {}));},
             [](PROfit::PROpeller &p, const py::buffer buf) {p.true_bin_indices = buffer_to_vector<int>(buf);}, 
             py::return_value_policy::reference_internal);

    // SystStruct
    py::class_<PROfit::SystStruct>(m, "SystStruct")
        .def(py::init<const std::string&, const int>())
        .def(py::init<const PROfit::SystStruct&>(), py::return_value_policy::copy)
        .def(py::init(&init_SystStruct_np))
        .def("GetSysName", &PROfit::SystStruct::GetSysName)
        .def("GetNUniverse", &PROfit::SystStruct::GetNUniverse)
        .def("HasWeightFormula", &PROfit::SystStruct::HasWeightFormula)
        .def("GetWeightFormula", &PROfit::SystStruct::GetWeightFormula)
        .def("SetWeightFormula", &PROfit::SystStruct::SetWeightFormula)
        .def("SetMode", &PROfit::SystStruct::SetMode)
        .def("CreateSpecs", &PROfit::SystStruct::CreateSpecs)
        .def("SanityCheck", &PROfit::SystStruct::SanityCheck)
        .def("FillCV", py::vectorize(&PROfit::SystStruct::FillCV))
        .def("FillUniverse", py::vectorize(&PROfit::SystStruct::FillUniverse))
        // CV and Variation must return the shared pointer since we manage ProSpec with the shared_ptr class
        .def("CV", [](PROfit::SystStruct& s) {return s.p_cv;})
        .def("Variation", [](PROfit::SystStruct& s, int universe) {return s.p_multi_spec.at(universe);})
        .def_readonly("systname",  &PROfit::SystStruct::systname)
        .def_readonly("n_univ",  &PROfit::SystStruct::n_univ)
        .def_readonly("mode",  &PROfit::SystStruct::mode)
        .def_readonly("weight_formula",  &PROfit::SystStruct::weight_formula)
        .def_readwrite("knobval", &PROfit::SystStruct::knobval)
        .def_readwrite("knob_index", &PROfit::SystStruct::knob_index)
        .def_readwrite("p_cv", &PROfit::SystStruct::p_cv)
        .def_readwrite("p_multi_spec", &PROfit::SystStruct::p_multi_spec)
        .def_readonly("index",  &PROfit::SystStruct::index);

    // PROsyst
    py::class_<PROfit::PROsyst>(m, "PROsyst")
        .def(py::init<>())
        .def(py::init<PROfit::PROsyst>())
        // PROsyst takes vector of systematics by reference. However, when passing a list from python,
        // There is no way to get around constructing a vector from the list (and copying SystStructs).
        // To make this explicit, we take the vector by value, and pass it to the class by reference.
        .def(py::init([](std::vector<PROfit::SystStruct> s) {return PROfit::PROsyst(s);}));
        // .def(py::init<const std::vector<PROfit::SystStruct>&>());

    // PROspec
    py::class_<PROfit::PROspec, std::shared_ptr<PROfit::PROspec>>(m, "PROspec")
        .def(py::init<>())
        .def(py::init<size_t>())
        .def(py::init<const PROfit::PROspec&>())
        .def("Spec", &PROfit::PROspec::Spec, py::return_value_policy::reference_internal) 
        .def("Error", &PROfit::PROspec::Error, py::return_value_policy::reference_internal);

    // PROsc
    py::class_<PROfit::PROsc>(m, "PROsc")
        .def(py::init<const PROfit::PROpeller&>())
        .def(py::init<PROfit::PROsc>());

    // instantiate numpy
    _import_array();
}
