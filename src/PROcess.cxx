#include "PROcess.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"
#include "PROtocall.h"

#include <Eigen/Eigen>

#include <random>
#include <string>
#include <vector>

namespace PROfit {
    PROspec FillCVSpectrum(const PROconfig &inconfig, const PROpeller &inprop, bool binned){

        PROspec myspectrum(inconfig.m_num_bins_total);

        if(binned) {
            for(size_t i = 0; i < inprop.hist.rows(); ++i) {
                float le = inprop.histLE[i];
                float systw = 1;
                for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                    myspectrum.Fill(k, systw * inprop.hist(i, k));
                }
            }
        } else {
            for(size_t i = 0; i<inprop.truth.size(); ++i){

                float add_w = inprop.added_weights[i]; 
                myspectrum.Fill(inprop.bin_indices[i], add_w);
            }
        }
        return myspectrum;

    }

    PROspec FillRecoSpectra(const PROconfig &inconfig, const PROpeller &inprop, const PROsc *inosc, const std::vector<float> &physparams, bool binned){

        PROspec myspectrum(inconfig.m_num_bins_total);

        if(binned) {
            for(long int i = 0; i < inprop.hist.rows(); ++i) {
                float le = inprop.histLE[i];
                if(physparams.size() != 0) {
                    for(size_t j = 0; j < inosc->model_functions.size(); ++j) {
                        float oscw = GetOscWeight(j, le, *inosc, physparams);
                        for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                            myspectrum.Fill(k, oscw * inosc->hists[j](i, k));
                        }
                    }
                } else {
                    for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                        myspectrum.Fill(k, inprop.hist(i, k));
                    }
                }
            }
        } else {
            for(size_t i = 0; i<inprop.truth.size(); ++i){
                float oscw  = physparams.size() != 0 ? GetOscWeight(i, inprop, *inosc, physparams) : 1;
                float add_w = inprop.added_weights[i]; 
                myspectrum.Fill(inprop.bin_indices[i], oscw * add_w);
            }
        }
        return myspectrum;

    }

    PROspec FillRecoSpectra(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst, const PROsc *inosc, const std::vector<float> &inshifts, const std::vector<float> &physparams, bool binned){

        PROspec myspectrum(inconfig.m_num_bins_total);

        if(binned) {
            for(long int i = 0; i < inprop.hist.rows(); ++i) {
                float le = inprop.histLE[i];
                float systw = 1;
                for(size_t j = 0; j < inshifts.size(); ++j) {
                    systw *= insyst.GetSplineShift(j, inshifts[j], i);
                }
                if(physparams.size() != 0) {
                    for(size_t j = 0; j < inosc->model_functions.size(); ++j) {
                        float oscw = GetOscWeight(j, le, *inosc, physparams);
                        for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                            myspectrum.Fill(k, systw * oscw * inosc->hists[j](i, k));
                        }
                    }
                } else {
                    for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                        myspectrum.Fill(k, systw * inprop.hist(i, k));
                    }
                }
            }
        } else {
            for(size_t i = 0; i<inprop.truth.size(); ++i){

                float oscw  = physparams.size() != 0 ? GetOscWeight(i, inprop, *inosc, physparams) : 1;
                float add_w = inprop.added_weights[i]; 
                const int true_bin = inprop.true_bin_indices[i]; 
                
                float systw = 1;
                for(size_t j = 0; j < inshifts.size(); ++j) {
                    systw *= insyst.GetSplineShift(j, inshifts[j], true_bin);
                }

                float finalw = oscw * systw * add_w;

                myspectrum.Fill(inprop.bin_indices[i], finalw);
            }
        }
        return myspectrum;

    }

  PROspec FillWeightedSpectrumFromHist(const PROconfig &inconfig, const PROpeller &inprop, const PROsc *inosc, std::vector<TH2D*> inweighthists, std::vector<float> &physparams, bool binned){
    PROspec myspectrum(inconfig.m_num_bins_total);

    if (binned) {
      log<LOG_WARNING>(L"%1% || WARNING: Binned fit is requested but not supported for histogram reweights. Returning empty spectrum so things will likely fail ") % __func__ ;
      // ETW Can't handle binned for now - will just return an empty spec
    }
    else {
      for(size_t i = 0; i<inprop.truth.size(); ++i){
	
	float oscw  = physparams.size() != 0 ? GetOscWeight(i, inprop, *inosc, physparams) : 1;
	float add_w = inprop.added_weights[i];
	float hist_w = 1.0 ;

	//Figure out what subchannel the event is in
	size_t subchan = inconfig.GetSubchannelIndexFromGlobalTrueBin(inprop.true_bin_indices[i]);
	std::string name = inconfig.m_fullnames[subchan];
	  
	//Put name for ICARUS study here. How to handle more generically?
	if (name == "nu_ICARUS_numu_numucc") {
	  double pmom = static_cast<double>(inprop.pmom[i]);
	  double pcosth = static_cast<double>(inprop.pcosth[i]);

	  for (size_t j = 0; j<inweighthists.size(); ++j){
	    TH2D h = *inweighthists[j];
	    int bin = h.FindBin(pmom,pcosth);
	    hist_w *= h.GetBinContent(bin);
	  }
	}

	float finalw = oscw * add_w * hist_w;
	log<LOG_DEBUG>(L"%1% || name %2% oscw %3% addw %4% histw %5%" ) % __func__ % name.c_str() %  oscw % add_w % hist_w;
	myspectrum.Fill(inprop.bin_indices[i], finalw);
      }
    }
    return myspectrum;
  }

    PROspec FillRecoSpectra(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst, const PROsc *inosc, const std::map<std::string, float> &inshifts, const std::vector<float> &physparams, bool binned){

        PROspec myspectrum(inconfig.m_num_bins_total);

        if(binned) {
            for(long int i = 0; i < inprop.hist.rows(); ++i) {
                float le = inprop.histLE[i];
                float systw = 1;
                for(const auto &[name, shift]: inshifts) {
                    systw *= insyst.GetSplineShift(name, shift, i);
                }
                if(physparams.size() != 0) {
                    for(size_t j = 0; j < inosc->model_functions.size(); ++j) {
                        float oscw = GetOscWeight(j, le, *inosc, physparams);
                        for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                            myspectrum.Fill(k, systw * oscw * inosc->hists[j](i, k));
                        }
                    }
                } else {
                    for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                        myspectrum.Fill(k, systw * inprop.hist(i, k));
                    }
                }
            }
        } else {
            for(size_t i = 0; i<inprop.truth.size(); ++i){

                float oscw  = physparams.size() != 0 ? GetOscWeight(i, inprop, *inosc, physparams) : 1;
                float add_w = inprop.added_weights[i]; 
                const int true_bin = inprop.true_bin_indices[i]; 
                
                float systw = 1;
                for(const auto &[name, shift]: inshifts) {
                    systw *= insyst.GetSplineShift(name, shift, true_bin);
                }

                float finalw = oscw * systw * add_w;

                myspectrum.Fill(inprop.bin_indices[i], finalw);
            }
        }
        return myspectrum;

    }

    PROspec FillRecoSpectra(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst, const std::map<std::string, float> &inshifts, bool binned) {

        PROspec myspectrum(inconfig.m_num_bins_total);

        if(binned) {
            for(long int i = 0; i < inprop.hist.rows(); ++i) {
                float systw = 1;
                for(const auto &[name, shift]: inshifts) {
                    systw *= insyst.GetSplineShift(name, shift, i);
                }
                for(size_t k = 0; k < myspectrum.GetNbins(); ++k) {
                    myspectrum.Fill(k, systw * inprop.hist(i, k));
                }
            }
        } else {
            for(size_t i = 0; i<inprop.truth.size(); ++i){
                float add_w = inprop.added_weights[i]; 
                const int true_bin = inprop.true_bin_indices[i]; 
                float systw = 1;
                for(const auto &[name, shift]: inshifts) {
                    systw *= insyst.GetSplineShift(name, shift, true_bin);
                }
                float finalw = systw * add_w;
                myspectrum.Fill(inprop.bin_indices[i], finalw);
            }
        }
        return myspectrum;
    }


    PROspec FillSystRandomThrow(const PROconfig &inconfig, const PROpeller &inprop, const PROsyst &insyst) {
        Eigen::VectorXf spec = Eigen::VectorXf::Constant(inconfig.m_num_bins_total, 0);
        Eigen::VectorXf cvspec = Eigen::VectorXf::Constant(inconfig.m_num_bins_total, 0);

        // TODO: We should think about centralizing rng in a thread-safe/thread-aware way
        static std::random_device rd{};
        static std::mt19937 rng{rd()};
        std::normal_distribution<float> d;
        std::vector<float> throws;
        //Eigen::VectorXf throwC = Eigen::VectorXf::Constant(inconfig.m_num_bins_total, 0);
        Eigen::VectorXf throwC = Eigen::VectorXf::Constant(inconfig.m_num_bins_total_collapsed, 0);
        for(size_t i = 0; i < insyst.GetNSplines(); i++)
            throws.push_back(d(rng));
        for(size_t i = 0; i < inconfig.m_num_bins_total_collapsed; i++)
            throwC(i) = d(rng);


        for(long int i = 0; i < inprop.hist.rows(); ++i) {
            float systw = 1;
            for(size_t j = 0; j < throws.size(); ++j) {
                systw *= insyst.GetSplineShift(j, throws[j], i);
            }
            for(size_t k = 0; k < inconfig.m_num_bins_total; ++k) {
                spec(k) += systw * inprop.hist(i, k);
                cvspec(k) += inprop.hist(i, k);
            }
        }

        if(insyst.GetNCovar() == 0) {
            Eigen::VectorXf final_spec = CollapseMatrix(inconfig, spec);
            return PROspec(final_spec, final_spec.array().sqrt());
        }

        // TODO: We probably just want to do this once and save it somewhere.
        // But where? PROpeller doesn't know about systs and PROsyst doesn't
        // know about cvspec.
        Eigen::MatrixXf diag = cvspec.asDiagonal();
        Eigen::MatrixXf full_cov = diag * insyst.fractional_covariance * diag;
        Eigen::LLT<Eigen::MatrixXf> llt(CollapseMatrix(inconfig, full_cov));

        Eigen::VectorXf final_spec = CollapseMatrix(inconfig, spec) + llt.matrixL() * throwC;
        
        return PROspec(final_spec, final_spec.array().sqrt());
    }

};
