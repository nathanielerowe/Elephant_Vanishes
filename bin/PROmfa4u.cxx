#include "PROconfig.h"
#include "PROspec.h"

#include "CLI11.h"
#include "LBFGSB.h"

#include <diy/master.hpp>
#include <diy/reduce.hpp>
#include <diy/partners/merge.hpp>
#include <diy/decomposition.hpp>
#include <diy/assigner.hpp>
#include <diy/mpi.hpp>
#include <diy/serialization.hpp>
#include <diy/partners/broadcast.hpp>
#include <diy/reduce-operations.hpp>
#include <diy/io/block.hpp>

#include <mfa/mfa.hpp>
#include <mfa/block_base.hpp>

#include <Eigen/Dense>
#include <Eigen/Eigen>
#include <Eigen/Dense>
#include <Eigen/SVD>
#include <Eigen/Core>
#include <unsupported/Eigen/CXX11/Tensor>
#include <unsupported/Eigen/NumericalDiff>

#define FMT_HEADER_ONLY
//#include <fmt/format.h>


using namespace PROfit;
log_level_t GLOBAL_LEVEL = LOG_DEBUG;

class ChiTest
{
    private:
        int n;
    public:
        ChiTest(int n_) : n(n_) {}
        float operator()(const Eigen::VectorXf &x, Eigen::VectorXf &grad)
        {
            float fx = 0.0;
            for(int i = 0; i < n; i += 2)
            {
                float t1 = 1.0 - x[i];
                float t2 = 10 * (x[i + 1] - x[i] * x[i]);
                grad[i + 1] = 20 * t2;
                grad[i]     = -2.0 * (x[i] * grad[i + 1] + t1);
                fx += t1 * t1 + t2 * t2;
            }
            return fx;
        }
};

// arguments to block foreach functions
struct DomainArgs
{
    DomainArgs(int dom_dim, int nvars) 
    {
        tot_ndom_pts = 0;
        starts.resize(dom_dim);
        ndom_pts.resize(dom_dim);
        full_dom_pts.resize(dom_dim);
        min.resize(dom_dim);
        max.resize(dom_dim);
        s.resize(nvars);
        f.resize(nvars);
        for (auto i = 0; i < nvars; i++)
        {
            s[i] = 1.0;
            f[i] = 1.0;
        }
        r = 0;
        t = 0;
        n = 0;
        multiblock = false;
        structured = true;   // Assume structured input by default
        rand_seed  = -1;
    }
    size_t              tot_ndom_pts;
    vector<int>         starts;                     // starting offsets of ndom_pts (optional, usually assumed 0)
    vector<int>         ndom_pts;                   // number of points in domain (possibly a subset of full domain)
    vector<int>         full_dom_pts;               // number of points in full domain in case a subset is taken
    vector<float>      min;                        // minimum corner of domain
    vector<float>      max;                        // maximum corner of domain
    vector<float>      s;                          // scaling factor for each variable or any other usage
    float              r;                          // x-y rotation of domain or any other usage
    vector<float>      f;                          // frequency multiplier for each variable or any other usage
    float              t;                          // waviness of domain edges or any other usage
    float              n;                          // noise factor [0.0 - 1.0]
    string              infile;                     // input filename
    bool                multiblock;                 // multiblock domain, get bounds from block
    bool                structured;                 // input data lies on unstructured grid
    int                 rand_seed;                  // seed for generating random data. -1: no randomization, 0: choose seed at random
};


// block
template <typename T>
struct Block : public BlockBase<T>
{
    using Base = BlockBase<T>;
    using Base::dom_dim;
    using Base::pt_dim;
    using Base::core_mins;
    using Base::core_maxs;
    using Base::bounds_mins;
    using Base::bounds_maxs;
    using Base::overlaps;
    using Base::input;

    static
        void* create()              { return mfa::create<Block>(); }

    static
        void destroy(void* b)       { mfa::destroy<Block>(b); }

    static
        void add(                                   // add the block to the decomposition
                int                 gid,                // block global id
                const Bounds<T>&    core,               // block bounds without any ghost added
                const Bounds<T>&    bounds,             // block bounds including any ghost region added
                const Bounds<T>&    domain,             // global data bounds
                const RCLink<T>&    link,               // neighborhood
                diy::Master&        master,             // diy master
                int                 dom_dim,            // domain dimensionality
                int                 pt_dim,             // point dimensionality
                T                   ghost_factor = 0.0) // amount of ghost zone overlap as a factor of block size (0.0 - 1.0)
        {
            mfa::add<Block, T>(gid, core, bounds, domain, link, master, dom_dim, pt_dim, ghost_factor);
        }

    static
        void save(const void* b_, diy::BinaryBuffer& bb)    { mfa::save<Block, T>(b_, bb); }
    static
        void load(void* b_, diy::BinaryBuffer& bb)          { mfa::load<Block, T>(b_, bb); }


    // read a floating point 2d scalar dataset from HDF5
    // reads masses for geometry dimension 0 from same HDF5 file
    // assigns integer values for the geometry dimension 1 from 0 to n_pts - 1
    // f = (mass, y, value)
    //template <typename V>               // type of science value being read
    void read_3d_data(
            const       diy::Master::ProxyWithLink& cp,
            //MFAInfo&    mfa_info,
            DomainArgs& args,
            Eigen::Tensor<float, 4> vals,
            bool  rescale)            // rescale science values
    {
        //std::cout << "$$$$ dom_dim: " << a->dom_dim << std::endl; 
        DomainArgs* a = &args;

        const int nvars       = 0;//mfa_info.nvars();
        const VectorXi mdims  = 0;//mfa_info.model_dims();

        // Resize the vectors that record error metrics
        this->max_errs.resize(nvars);
        this->sum_sq_errs.resize(nvars);

        // Set points per direction and compute total points in domain
        VectorXi ndom_pts(dom_dim);
        for (int i = 0; i < dom_dim; i++)
            ndom_pts(i)     =  a->ndom_pts[i];

        // Create input data set and add to block
        //std::cout << "dom_dim, mdims, ndom_pts.prod(), ndom_pts: " << dom_dim << ", " << mdims << ", " << ndom_pts.prod() << ", " << ndom_pts << std::endl;
        input = new mfa::PointSet<T>(dom_dim, mdims, ndom_pts.prod(), ndom_pts);
        assert(vals(0) == ndom_pts(0));
        assert(vals(1) == ndom_pts(1));
        assert(vals(2) == ndom_pts(2));
        // set geometry values
        int n = 0;
        int pd = 0;//mfa_info.pt_dim();
        for (size_t k = 0; k < (size_t)(ndom_pts(2)); k++) {
            for (size_t j = 0; j < (size_t)(ndom_pts(1)); j++) {
                for (size_t i = 0; i < (size_t)(ndom_pts(0)); i++) {
                    input->domain(n, 0) = i;
                    input->domain(n, 1) = j;
                    input->domain(n, 2) = k;
                    for (int m = 0; m < pd-dom_dim; m++) {
                        //std::cout << "index n, 3+m, i, j, k, vals: " << n << ", " << 3+m << ", " << i << ", " << j << ", " << k << ", ";
                        input->domain(n, dom_dim+m) = vals(m, i, j, k);
                        //std::cout << vals(m, i, j, k) << std::endl;

                    }	
                    n++;
                }
            }
        }

        // Init params from input data (must fill input->domain first)
        input->init_params();   

        // Construct the MFA object
        //this->setup_MFA(cp, mfa_info);


        // find extent of masses, values, and science variables (bins)
        // NOTE: the difference between core_mins/maxs and bounds_mins/maxs only matters
        //       if DIY is used to partition the domain space, but we set them appropriately anyway
        bounds_mins.resize(pt_dim);
        bounds_maxs.resize(pt_dim);
        core_mins.resize(dom_dim);
        core_maxs.resize(dom_dim);
        bounds_mins = input->domain.colwise().minCoeff();
        bounds_maxs = input->domain.colwise().maxCoeff();
        core_mins = bounds_mins.head(dom_dim);
        core_maxs = bounds_maxs.head(dom_dim);

        //std::cout << "tot_ndom_pt, input->domain(tot_ndom_pts - 1, 1), this->dom_dim = " << tot_ndom_pts << ", " << input->domain(tot_ndom_pts - 1, 1) << ", " << this->dom_dim << std::endl;

        // debug
        //cerr << "domain extent:\n min\n" << this->bounds_mins << "\nmax\n" << this->bounds_maxs << endl;
    }

};

int main(int argc, char* argv[])
{

    CLI::App app{"Test for PROfit"}; 

    // Define options
    std::string xmlname = "NULL.xml"; 
    int maxevents = 100;

    //floats
    app.add_option("-x,--xml", xmlname, "Input PROfit XML config.");
    app.add_option("-m,--max", maxevents, "Max number of events to run over.");
    app.add_option("-v,--verbosity", GLOBAL_LEVEL, "Verbosity Level [1-4].");

    CLI11_PARSE(app, argc, argv);


    diy::mpi::environment env(argc, argv);
    diy::mpi::communicator world;

    std::cout<<"World Rank: "<<world.rank()<<std::endl;
    size_t blocks = world.size();
    if (world.rank()==0) std::cout<<"We have blocks: "<< blocks<<std::endl;

    int mfadim = 5; //Dimension of the MFA model
    Bounds<float> fc_domain(3);
    fc_domain.min[0] = 0.;
    fc_domain.max[0] = blocks-1;
    fc_domain.min[1] = 0.;
    fc_domain.max[1] = blocks-1;
    fc_domain.min[2] = 0.;
    fc_domain.max[2] = blocks-1;

    diy::FileStorage               storage("./DIY.XXXXXX");
    diy::RoundRobinAssigner        fc_assigner(world.size(), blocks);
    diy::RegularDecomposer<Bounds<float>> fc_decomposer(mfadim, fc_domain, blocks);
    diy::RegularBroadcastPartners  fc_comm(    fc_decomposer, mfadim, true);
    diy::RegularMergePartners      fc_partners(fc_decomposer, mfadim, true);
    diy::Master                    fc_master(world, 1, -1, &Block<float>::create, &Block<float>::destroy, &storage, &Block<float>::save, &Block<float>::load);
    diy::ContiguousAssigner   assigner(world.size(), blocks);



    PROconfig myConf(xmlname);

    //PROspec mySpec(myConf);
    //TH1D hmm = mySpec.toTH1D(myConf);


    LBFGSpp::LBFGSBParam<float> param;  
    param.epsilon = 1e-6;
    param.max_iterations = 100;
    LBFGSpp::LBFGSBSolver<float> solver(param); 

    int n=78;
    ChiTest fun(n);

    // Bounds
    Eigen::VectorXf lb = Eigen::VectorXf::Constant(n, 0.0);
    Eigen::VectorXf ub = Eigen::VectorXf::Constant(n, std::numeric_limits<float>::infinity());

    // Initial guess
    Eigen::VectorXf x = Eigen::VectorXf::Constant(n, 2.0);


    // x will be overwritten to be the best point found
    float fx;
    int niter = solver.minimize(fun, x, fx, lb, ub);


    std::cout << niter << " iterations" << std::endl;
    std::cout << "x = \n" << x.transpose() << std::endl;
    std::cout << "f(x) = " << fx << std::endl;

    return 0;
}

