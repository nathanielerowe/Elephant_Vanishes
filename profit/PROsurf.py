import argparse
import profit
import pandas as pd
import uproot

def main(xml=""):
    c = profit.PROconfig(xml)
    print("DONE PARSING")

    # process files
    num_files = c.m_num_mcgen_files

    for f, t, bs in zip(c.m_mcgen_file_name, c.m_mcgen_tree_name, c.m_branch_variables):
        if f.endswith(".root"):
             process_root_file(f, t, bs, c)

def process_root_file(fname, tree, branches, c):
    tf = uproot.open(fname)
    ttree = tf[tree]

    friends = []
    for friend, friend_f in zip(c.m_mcgen_file_friend_treename_map.get(fname, []), c.m_mcgen_file_friend_map.get(fname, [])):
        friends.append(uproot.open(friend_f)[friend])

    # TODO: include other friends
    for ib in range(len(branches)):
        branch = branches[ib]
        branches[ib].branch_formula = "branch_form_%s_%i" % (fname, ib), branch.name, ttree


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xml", type=str, help="XML file name.", default="")
    return parser.parse_args().__dict__

if __name__ == "__main__":
    args = parse_args()
    main(**args)
