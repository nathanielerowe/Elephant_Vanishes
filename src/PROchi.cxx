#include "PROchi.h"
#include "Eigen/src/Core/Matrix.h"
#include "PROlog.h"
using namespace PROfit;


PROchi::PROchi(const std::string tag, const PROconfig *conin, const PROpeller *pin, const PROsyst *systin, const PROsc *oscin, const PROspec &datain, int nparams) : model_tag(tag), config(conin), peller(pin), syst(systin), osc(oscin), data(datain), nparams(nparams) {
    last_value = 0.0; last_param = Eigen::VectorXd::Zero(nparams); 
}

float PROchi::operator()(const Eigen::VectorXd &param, Eigen::VectorXd &gradient){

    //std::vector<float> shifts = param(Eigen::SeqN(2,1)).array();
    //std::vector<float> fitparams = param(Eigen::SeqN(0,2)).array();

    // Get Spectra from FillRecoSpectra
    Eigen::VectorXd subvector1 = param.segment(0, 2);
    std::vector<float> fitparams(subvector1.data(), subvector1.data() + subvector1.size());
    Eigen::VectorXd subvector2 = param.segment(2,nparams - 2);
    std::vector<float> shifts(subvector2.data(), subvector2.data() + subvector2.size());

    log<LOG_DEBUG>(L"%1% || Shifts size is %2%") % __func__ % shifts.size();

    PROspec result = FillRecoSpectra(*config, *peller, *syst, *osc, shifts, fitparams);

    std::cout<<"Spec "<< result.Spec()<<" .. "<<std::endl;
    result.Print();

    // Calcuate Full Covariance matrix
    Eigen::MatrixXd diag = result.Spec().array().matrix().asDiagonal(); 
    Eigen::MatrixXd full_covariance =  diag*(syst->fractional_covariance)*diag;
    std::cout<<"Full: "<<full_covariance.size()<<std::endl;
    std::cout<<full_covariance<<std::endl;

    // Collapse Covariance and Spectra 
    Eigen::MatrixXd collapsed_full_covariance =  CollapseMatrix(*config,full_covariance);  
    std::cout<<"cFull: "<<collapsed_full_covariance.size()<<std::endl;
    std::cout<<collapsed_full_covariance<<std::endl;

    Eigen::MatrixXd stat_covariance = data.Spec().array().matrix().asDiagonal();
    Eigen::MatrixXd collapsed_stat_covariance = CollapseMatrix(*config, stat_covariance); 
    std::cout<<"cStat: "<<collapsed_stat_covariance.size()<<std::endl;
    std::cout<<collapsed_stat_covariance<<std::endl;

    log<LOG_DEBUG>(L"%1% || HERE") % __func__ ;

    // Invert Collaped Matrix Matrix 
    Eigen::MatrixXd inverted_collapsed_full_covariance = (collapsed_full_covariance+collapsed_stat_covariance).inverse();

    std::cout<<"shape: "<<inverted_collapsed_full_covariance.size()<<std::endl;
    std::cout<<inverted_collapsed_full_covariance<<std::endl;

    for (int i = 0; i < nparams; i++) {
        if(param(i) == last_param(i)) {
            gradient(i) = 1;
            continue;
        }
        Eigen::VectorXd tmpParams = last_param;
        tmpParams(i) = param(i);
        //Eigen::VectorXd subvector2 = tmpParams;
        Eigen::VectorXd subvector1 = tmpParams.segment(0, 2);
        std::vector<float> fitparams(subvector1.data(), subvector1.data() + subvector1.size());
        Eigen::VectorXd subvector2 = tmpParams.segment(2,nparams - 2);
        std::vector<float> shifts(subvector2.data(), subvector2.data() + subvector2.size());
        PROspec result = FillRecoSpectra(*config, *peller, *syst, *osc, shifts, fitparams);
        // Calcuate Full Covariance matrix
        Eigen::MatrixXd diag = result.Spec().array().matrix().asDiagonal(); 
        Eigen::MatrixXd full_covariance =  diag*(syst->fractional_covariance)*diag;
        // Collapse Covariance and Spectra 
        Eigen::MatrixXd collapsed_full_covariance =  CollapseMatrix(*config,full_covariance);  
        Eigen::MatrixXd stat_covariance = data.Spec().array().matrix().asDiagonal();
        Eigen::MatrixXd collapsed_stat_covariance = CollapseMatrix(*config, stat_covariance); 
        // Invert Collaped Matrix Matrix 
        Eigen::MatrixXd inverted_collapsed_full_covariance = (collapsed_full_covariance+collapsed_stat_covariance).inverse();
        // Calculate Chi^2  value
        Eigen::VectorXd delta  = result.Spec() - data.Spec(); 
        float pull = subvector2.array().square().sum(); 
        float value = (delta.transpose())*inverted_collapsed_full_covariance*(delta) + pull;
        gradient(i) = (value-last_value)/(param(i) - last_param(i));
    }

    std::cout<<"Grad: "<<gradient<<std::endl;

    // Calculate Chi^2  value
    Eigen::VectorXd delta  = result.Spec() - data.Spec(); 
    float pull = subvector2.array().square().sum(); 
    float value = (delta.transpose())*inverted_collapsed_full_covariance*(delta) + pull;

    //log<LOG_DEBUG>(L"%1% || value %2%, last_value %3%, pull") % __func__ % value  % last_value % pull;
    log<LOG_DEBUG>(L"%1% || FINISHED ITERATION got vals: %2% %3%") % __func__ % value % last_value ;

    //Update last param
    last_param = param;
    last_value = value;

    return value;
}


