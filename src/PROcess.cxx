#include "PROcess.h"
#include "PROlog.h"
#include "PROspec.h"
#include "PROsyst.h"

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
	double pmom = static_cast<double>(inprop.pmom[i]);
	double pcosth = static_cast<double>(inprop.pcosth[i]);
	double pdg = static_cast<double>(inprop.pdg[i]);

	for (size_t j = 0; j<inweighthists.size(); ++j){
	  TH2D h = *inweighthists[j];
	  int bin = h.FindBin(pmom,pcosth);
	  hist_w *= h.GetBinContent(bin);
	}
	
	float finalw = oscw * add_w * hist_w;
	log<LOG_DEBUG>(L"%1% || PDG: %2% Mom: %3% Costh %3% Hist_w %5% ") % __func__ % pdg % pmom % pcosth % hist_w;	
	log<LOG_DEBUG>(L"%1% || Oscw: %2% Add_w %3% Hist_w %4% Finalw %5%") % __func__ % oscw % add_w % hist_w % finalw;	
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

  float GetOscWeight(int rule, float le, const PROsc &inosc, const std::vector<float> &inphysparams) {
      
        // The model functions take L and E separately, so give E as 1 and L as L/E
        return inosc.model_functions[rule](std::pow(10, inphysparams[0]), std::pow(10, inphysparams[1]), 1, le);
    }

  float GetOscWeight(int ev_idx, const PROpeller &inprop, const PROsc &inosc, const std::vector<float> &inphysparams) {

        //get subchannel here from pdg. this will be added when we agree on convention.
        //for now everything is numu disappearance 3+1. 
        // inphysparams[0] is log(delta-msq)
        // inphysparams[1] is sinsq2thmumu

//        float prob2 = inosc.Pmumu(std::pow(10, inphysparams[0]), inphysparams[1], inprop.truth[ev_idx], inprop.baseline[ev_idx]);
        float prob = inosc.model_functions[inprop.model_rule[ev_idx]](std::pow(10, inphysparams[0]), std::pow(10, inphysparams[1]), inprop.truth[ev_idx], inprop.baseline[ev_idx] );
        return prob;
    }

};
