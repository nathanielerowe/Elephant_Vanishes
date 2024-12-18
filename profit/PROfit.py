import argparse
import profit
import pandas as pd

def main(xml=""):
    c = profit.PROconfig(xml)
    print("DONE PARSING")

    # process files
    num_files = c.m_num_mcgen_files

    for f in c.m_mcgen_file_name:
        df = pd.read_hdf(f, keyname="df")
        print(f)

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xml", type=str, help="XML file name.", default="")
    return parser.parse_args().__dict__

if __name__ == "__main__":
    args = parse_args()
    main(**args)
