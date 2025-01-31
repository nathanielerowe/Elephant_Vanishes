![Alt text](/other/profit.png "Minimizing PROfit")


## Installing through Python Interface

You should start from an environment with ROOT, HDF5, and BOOST already installed or setup. These can be installed on linux with apt-get, or on mac with homebrew. On the FNAL servers, you can get these setup by running:
```
source /cvmfs/larsoft.opensciencegrid.org/setup_larsoft.sh
setup root v6_28_12 -q e26:p3915:prof
setup cmake v3_27_4
setup hdf5 v1_12_2b -q e26:prof
setup boost v1_82_0 -q e26:prof
```

With these packages avialble, setup a python environment:
```
python -m venv env
. env/bin/activate
```
Install the build dependencies for profit:
```
pip install --upgrade pip
pip install wheel setuptools pybind11 numpy==2.0.2
```
Install profit:
```
pip install git+https://github.com/gputnam/Elephant_Vanishes
```

Then you're done! You can now ``import profit`` in a python shell, run the python executables (``PROsurf.py``, e.g.), or run the PRO* binary executables with the ``PRO`` helper (``PRO PROsurf``, e.g.).

### Installing for Development

If you want to develop PROfit, then it is useful to clone down a local copy of the repository and then install it. Start with the same steps to initialize your python virtual environment:
```
python -m venv env
. env/bin/activate
pip install --upgrade pip
pip install wheel setuptools pybind11 numpy===2.0.2
```

Then pull down the profit repository:
```
git clone https://github.com/gputnam/Elephant_Vanishes.git
```
Then enter the repsoitory and install it:
```
cd Elephant_Vanishes
pip install .
```
If you make some changes, you can update the installation by re-running ``pip install .``

## Usage

### Basic Description of main classes and use cases

- **PROlog** : Fairly simple logging/verbosity wrapper class. This is the only way PROfit should print. Usage as:
    -       log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;
    - LOG CRITICAL = 0
    - LOG ERROR = 1
    - LOG WARNING = 2
    - LOG INFO = 3 :  This should be level of standard physics output that someone would plausibly use, but not create giant log files!
    - LOG DEBUG = 4 
- **PROconfig** : Primary bookkeeping and XML loading class. Stores ALL information of mode-detector-channel-subchannel information. Must be created once and only once per PROfit executable and passed by reference to more complex functions when it is needed.
    - Contains helper functions channel<-> subchannel mapping..etc..
    - As well as functions for how to collapse subchannels to channels
- **PROspec**  : Class for storing final spectra and errors. Barebones class as all binning,collapsing handeld by PRoconfig/
    - Essentially two Eigen::VectorXD ! Simple with helper functions
- **PROpeller** : The PROpeller, which moves the analysis forward. A class to keep all MC events for oscllation event-by-event.
    - Saved as a set of six std::vector floats or ints. 
    - Currently only truth (energy), reco (energy) baseline, pdg, addedweight and precalculated bin incicies indicating where the reco variable fits in a given PROconfig defined binning. 
- **PROtocall** : Currently a bit of a rogue set of functions that, given a PROconfig, dictate how to get binning and collapse Covariance matricies
- **PROcess** : The master weighting function that combines everything to give one final output spectra (PROspec) event-by-event. 
    - FillRecoSpectra (PROconfig,PROpeller,PROsyst,PROosc, shifts, phys param)
- **PROcreate** : Workhorse classes for loading the set of root ntuples and systematics into PROspec CV class and a vector of SystStructs for systematics.  
    - Needs a PROconfig to tell it exactly what to load and where (channel and subchannel info) to assign the variables. Also which systeatics to load.
    - The class **SystStruct** is also defined in here (for some reason), which stores the weights and spline information for systematics
- **PROsyst** : Class that groups all systematics (each with a SystStruct) and manages their formation and effect on PROspecs
    - i.e generates covariance matricies, fills splines, gets spline shifted spectrum.. etc..
    - Currently has Spline (CAFana style), Covariance (SBNfit style) and MFA (2D spline to-be-implemented) options
- **PROsc** : The model class (not the lack of second 'o'). Currently is a very simple hardcoded 3+1 SBL approximation.
    -- End goal to implement custom classes, NuSquids..etc..
- **PROchi** : Class that gathers the MC (PROpeller), Systematics (PROsyst) and model (PROsc) and forms a function calculating a chi^2 that can be minimized over


### Some Possible PRO class names abailable 
PROtest
PROfessional
PROspect
PROton
PROpane
PROcrastiation
PRO.file
PROtect
PROcreate
PROduce
PROgress
PROfane
PROfanity
PROlapse
PROtocoll


### Junk notes below here. Put somewhere better 
  rec.slc.truth.wgt..length = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 0, 0                                                         
  rec.slc.truth.wgt..totarraysize = 130                                                                
  rec.slc.truth.wgt.univ..length = 6,  6, 6, 6, 6, 6, 6, 6, 6, 6, 6,  6, 6, 6, 6, 6, 6, 6, 6, 6, rec.slc.truth.wgt.univ..length[rec.slc.truth.wgt..totarraysize
  rec.slc.truth.wgt.univ..totarraysize = 17559
  rec.slc.truth.wgt.univ = 1,  1, 1, 1, 1, 1,   1, 1, 1, 1, 1,   1, 1, 1, 1, 1, 1, 1, 0.861389, 1.15461   .. of length [rec.slc.truth.wgt.univ..totarraysize - 17559 Float] 
  rec.slc.truth.wgt.univ..idx = 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114  .. of lngth [rec.slc.truth.wgt..totarraysize INT]
  rec.slc.truth.wgt..idx = 0,            of length [ rec.slc..length ]

  (rec.slc[si].truth.wgt[wi].univ[take next N])    [NUM]
  rec.slc.truth.wgh.univ[rec.slc.truth.wgt..idx[SLC]+I]+6
   

