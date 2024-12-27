import argparse
import profit

# TODO: WHY DOES SYSTSTRUCT USE unique_ptr's TO CV + SYST UNIVERSES???

def main(xml="", inject=None):
    profit.Globals.GLOBAL_LEVEL = verbosity

    # Init the config object
    c = profit.PROconfig(xml)
    profit.PROlogINFO("DONE PARSING")

    # Load the SystStruct's and PROpeller
    # 
    # This can be used interchangably with  profit.PROcess_CAFana
    syst_structs, prop = profit.PROcess_dataframes(c)

    # Load systematic weights into PROsyst
    systs = profit.PROsyst(syst_structs)

    # Define oscillation model
    osc = profit.PROsc(prop)

    profit.PROlogINFO("DONE LOADING PRO CLASSES")

    # Injected model point
    return

def two_floats(s):
    return float(s.split(" ")[0]), float(s.split(" ")[1])

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xml", type=str, help="XML file name.", default="")
    parser.add_argument("-i", "--inject", type=two_floats, help="Physics parameters to inject as true signal.", default="0 0")
    return parser.parse_args().__dict__

if __name__ == "__main__":
    args = parse_args()
    main(**args)
