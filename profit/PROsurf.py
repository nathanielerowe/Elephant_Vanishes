import argparse
import profit

# TODO: WHY DOES SYSTSTRUCT USE unique_ptr's TO CV + SYST UNIVERSES???

def main(xml="", 
    inject=None, inject_systs={}, 
    eventbyevent=False, 
    systs_list=None, exclude_systs=None, 
    grid=(40,), linx=False, liny=False,
    xlo=1e-4, xhi=1, ylo=1e-2, yhi=1e2,
    statonly=False, outfile="", nthread=1,
    xlabel="", ylabel="", verbosity=0):

    # additional config
    logx = not linx
    logy = not liny
    savetoroot = outfile.endswith(".root")
    profit.Globals.GLOBAL_LEVEL = verbosity

    # Init the config object
    c = profit.PROconfig(xml)
    profit.PROlogINFO("DONE PARSING")

    # Load the SystStruct's and PROpeller
    # 
    # This can be used interchangably with profit.PROcess_CAFana
    syst_structs, prop = profit.PROcess_dataframes(c)

    # Load systematic weights into PROsyst
    systs = profit.PROsyst(syst_structs)

    # Define oscillation model
    osc = profit.PROsc(prop)

    profit.PROlogINFO("DONE LOADING PRO CLASSES")

    # Information on parameter injection
    profit.PROlogINFO("Injected point: sinsq2th = %f dmsq = %f" % (inject[0], inject[1]))
    for name, shift in inject_systs.items(): 
        profit.PROlogINFO("Injected syst: %s shifted by %f sigma" % (name, shift))

    # Asimov data
    data = None
    if inject[0] != 0 and inject[1] != 0:
        pparams = (np.log10(inject[0]), np.log10(inject[1]))
        data = profit.FillRecoSpectra(c, prop, systs, osc, inject_systs, pparams, not eventbyevent)
    elif len(inject_systs):
        data = profit.FillRecoSpectra(c, prop, systs, inject_systs, not eventbyevent)
    else:
        data = profit.FillCVSpectrum(c, prop, not eventbyevent)
    
    if systs_list is not None:
        systs = systs.subset(systs_list)
    elif exclude_systs is not None:
        systs = systs.excluding(exclude_systs)

    if len(grid) == 1:
        nbinsx = grid[0]
        nbinsy = grid[0]
    else:
        nbinsx = grid[0]
        nbinsy = grid[1]

    # Compute the chi2 surface
    surface = profit.PROsurf(nbinsx, profit.LogLin.LogAxis if logx else profit.LogLin.LinAxis, xlo, xhi,
                             nbinsy, profit.LogLin.LogAxis if logy else profit.LogLin.LinAxis, ylo, yhi)

    if statonly:
        surface.FillSurfaceStat(c, prop, osc, data, outfile if not savetoroot else "", not eventbyevent)
    else:
        surface.FillSurface(c, prop, systs, osc, data, outfile if not savetoroot else "", not eventbyevent, nthread)

    # Detect if we are saving to ROOT file
    rootfile = outfile if savetoroot else ""

    # Write output
    profit.savePROsurf(surface, logx, logy, xlabel, ylabel, rootfile, "PROfit_surface.pdf")

def tupling(*types):
    def get(s):
        dat = s.split(",")
        return [t(d) for t, d in zip(types, dat)]
    return get

class ParseDict(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        """Parses key-value pairs from a list of strings and stores them in a dictionary."""
        dict_value = {}
        for key, val in values:
            dict_value[key] = val
        setattr(namespace, self.dest, dict_value)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xml", type=str, help="XML file name.", default="")
    parser.add_argument("-i", "--inject", type=tupling(float, float), help="Physics parameters to inject as true signal.", default="0,0")
    parser.add_argument("--inject-systs", type=tupling(str, float), action=ParseDict, help="Systematic shifts to inject. Map of name and shift value in sigmas. Only spline systs are supported right now.", default={}, nargs="*")
    parser.add_argument("--eventbyevent", action="store_true", help="Do you want to weight event-by-event?")
    parser.add_argument("--systs-list", nargs="*", default=None, help="Override list of systematics to use (note: all systs must be in the xml).")
    parser.add_argument("--exclude-systs", nargs="*", default=None, help="List of systematics to exclude.")
    parser.add_argument("--grid", type=tupling(int, int), default=(40,))
    parser.add_argument("--linx", action="store_true", help="Specify if x-axis is logarithmic or linear (default log)")
    parser.add_argument("--liny", action="store_true", help="Specify if y-axis is logarithmic or linear (default log)")
    parser.add_argument("--xlo", type=float, default=1e-4, help="Lower limit for x-axis")
    parser.add_argument("--xhi", type=float, default=1, help="Upper limit for x-axis")
    parser.add_argument("--ylo", type=float, default=1e-2, help="Lower limit for y-axis")
    parser.add_argument("--yhi", type=float, default=1e2, help="Upper limit for y-axis")
    parser.add_argument("--statonly", action="store_true", help="Run a stats only surface instead of fitting systematics")
    parser.add_argument("-o", "--outfile", help="If you want chisq to be dumped to text file, provide name", default="")
    parser.add_argument("-t", "--nthread", type=int, default=1, help="Number of fits.")
    parser.add_argument("--xlabel", type=str, default="sin^{2}2#theta_{#mu#mu}", help="X-axis label")
    parser.add_argument("--ylabel", type=str, default="#Deltam^{2}_{41}", help="Y-axis label")
    parser.add_argument("-v", "--verbosity", type=int, default=profit.Globals.GLOBAL_LEVEL, help="Verbosity Level [1-4].")
    return parser.parse_args().__dict__

if __name__ == "__main__":
    args = parse_args()
    main(**args)
