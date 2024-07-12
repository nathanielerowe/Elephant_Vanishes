![Alt text](/other/profit.png "Minimizing PROfit")







### Basic Description of main classes and use cases

- **PROlog** : Fairly simple logging/verbosity wrapper class. This is the only way PROfit should print. Usage as:
    -       log<LOG_DEBUG>(L"%1% || Using a total of %2% individual files") % __func__  % num_files;
    - LOG_CRITICAL = 0
    - LOG_ERROR = 1
    - LOG_WARNING = 2
    - LOG_INFO = 3 :  This should be level of standard physics output that someone would plausibly use, but not create giant log files!
    - LOG_DEBUG = 4 
- **PROconfig** : Primary bookkeeping and XML loading class. Stores ALL information of mode_detector_channel_subchannel information. Must be created once and only once per PROfit executable and passed by reference to more complex functions when it is needed.
    - Contains helper functions channel<-> subchannel mapping..etc..
    - As well as functions for how to collapse subchannels->channels
- **PROspec**  : Class for storing final spectra and errors. Barebones class as all binning,collapsing handeld by PRoconfig/
    - Essentially two Eigen::VectorXD ! Simple with helper functions
- **PROpeller** : The PROpeller, which moves the analysis forward. A class to keep all MC events for oscllation event-by-event.
    - Saved as a set of six std::vector floats or ints. 
    - Currently only truth (energy), reco (energy) baseline, pdg, added_weight and precalculated bin incicies indicating where the reco variable fits in a given PROconfig defined binning. 
- **PROtocall** : Currently a bit of a rogue set of functions that, given a PROconfig, dictate how to get binning and collapse Covariance matricies
- **SystStruct** : A 
- **PROsyst** : Class that groups all systematics (each with a SystStruct) and manages their formation and effect on PROspecs
    - i.e generates covariance matricies, fills splines, gets spline shifted spectrum.. etc..
    - Currently 
- **PROchi** : Class that gathers the MC (PROpeller), Systematics (PROsyst) and model (PROsc) and forms a function calculating a chi^2 that can be minimized over




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
#PROtocoll

  rec.slc.truth.wgt..length = 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 130, 0, 0                                                         
  rec.slc.truth.wgt..totarraysize = 130                                                                
  rec.slc.truth.wgt.univ..length = 6,  6, 6, 6, 6, 6, 6, 6, 6, 6, 6,  6, 6, 6, 6, 6, 6, 6, 6, 6, rec.slc.truth.wgt.univ..length[rec.slc.truth.wgt..totarraysize
  rec.slc.truth.wgt.univ..totarraysize = 17559
  rec.slc.truth.wgt.univ = 1,  1, 1, 1, 1, 1,   1, 1, 1, 1, 1,   1, 1, 1, 1, 1, 1, 1, 0.861389, 1.15461   .. of length [rec.slc.truth.wgt.univ..totarraysize - 17559 Float] 
  rec.slc.truth.wgt.univ..idx = 0, 6, 12, 18, 24, 30, 36, 42, 48, 54, 60, 66, 72, 78, 84, 90, 96, 102, 108, 114  .. of lngth [rec.slc.truth.wgt..totarraysize INT]
  rec.slc.truth.wgt..idx = 0,            of length [ rec.slc..length ]

  (rec.slc[si].truth.wgt[wi].univ[take next N])    [NUM]
  rec.slc.truth.wgh.univ[rec.slc.truth.wgt..idx[SLC]+I]+6
   



TODO:
Check const for all functions that should have const!
Jenkins hash the co?

