#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <numpy/arrayobject.h>
#include <string>
#include <format>
#include <iostream>
#include "PROconfig.h"
#include "PROcreate.h"
#include "TTreeFormula.h"

namespace py = pybind11;

// expected by PROfit for printing stuff
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

template <typename T>
std::vector<T> buffer_to_vector(const py::buffer buf) {
    py::buffer_info info = buf.request();

    // Check if the buffer is the correct type and dimensions
    if (info.ndim > 1) {
        throw std::runtime_error("Buffer must be one-dimensional!");
    }

    // Special case: convert buffer of doubles to vector of floats
    if (info.format == py::format_descriptor<double>::format() &&
        py::format_descriptor<T>::format() == py::format_descriptor<float>::format()) {

      std::vector<float> vec(info.size);
      double* ptr = static_cast<double*>(info.ptr);

      for (ssize_t i = 0; i < info.size; ++i) {
          vec[i] = static_cast<float>(ptr[i]);
      }

      return vec;

    }
    if (info.format != py::format_descriptor<T>::format()) {
        std::stringstream err;
        err << "Incompatible buffer format! Expected (" << py::format_descriptor<T>::format() << "), received (" <<  info.format <<  ").";
        throw std::runtime_error(err.str());
    }

    std::vector<T> vec(info.size);
    T* ptr = static_cast<T*>(info.ptr);

    // Copy data from buffer to vector
    std::copy(ptr, ptr + info.size, vec.begin());

    return vec;
}

std::string ttreeformula_getter(const TTreeFormula *f) {
  if (!f) return "NULL";

  return std::string(f->PrintValue());
}

PROfit::SystStruct _init_SystStruct_np(
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

    py::class_<PROfit::BranchVariable, std::shared_ptr<PROfit::BranchVariable>>(m, "BranchVariable")
        .def(py::init<PROfit::BranchVariable>())
        .def(py::init<std::string, std::string, std::string>())
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

    py::class_<PROfit::PROconfig>(m, "PROconfig")
        .def(py::init<const std::string&>())
        .def(py::init<>())
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

    py::class_<PROfit::SystStruct>(m, "SystStruct")
        .def(py::init<const std::string&, const int>())
        .def(py::init(&_init_SystStruct_np));

    // m.def("init_PROpeller", &_init_PROpeller, "Initialize PROpeller object");

    // instantiate numpy
    _import_array();
}
