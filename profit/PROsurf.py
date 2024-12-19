import argparse
import profit
import numpy as np
import pandas as pd
import uproot
import awkward as ak

def main(xml=""):
    c = profit.PROconfig(xml)
    print("DONE PARSING")

    # process files
    num_files = c.m_num_mcgen_files

    for f, t, bs in zip(c.m_mcgen_file_name, c.m_mcgen_tree_name, c.m_branch_variables):
        if f.endswith(".root"):
             process_root_file(f, t, bs, c)

def make_array(ttree, branches):
    return ak.to_dataframe(ttree.arrays(branches, library="ak"))

def process_root_file(fname, tree, branches, c):
    tf = uproot.open(fname)
    ttree = tf[tree]

    friends = []
    friend_names = []
    for friend, friend_f in zip(c.m_mcgen_file_friend_treename_map.get(fname, []), c.m_mcgen_file_friend_map.get(fname, [])):
        # Remove duplicates (TODO: why are they there???)
        if friend in friend_names: continue

        friend_names.append(friend)
        friends.append(uproot.open(friend_f)[friend])

    # TODO: include other friends

    # get the event weights
    event_weights_pertree = [make_array(t, [k for k in t.keys() if k in c.m_mcgen_variation_allowlist and k not in c.m_mcgen_variation_denylist]) 
        for t in [ttree] + friends]

    # map_systematic_num_universe <- number of universes for each systematic
    map_systematic_num_universe = {}
    for ew in event_weights_pertree:
        if ew is not None:
            for col in ew.columns:
                map_systematic_num_universe[col] = ew.index.get_level_values(1).max()

    for ib in range(len(branches)):
        branch = branches[ib]
        branches[ib].branch_formula = "branch_form_%s_%i" % (fname, ib), branch.name, ttree

    # List of systematics
    syst_structs = []

    total_num_systematics = len(map_systematic_num_universe) 
    for syst, nuniv in map_systematic_num_universe.items():
        mode = c.m_mcgen_variation_type_map[syst]

        syst_structs.append(profit.SystStruct(syst, nuniv))
        syst_structs[-1].SetWeightFormula("1")
        syst_structs[-1].SetMode(mode)

        if mode == "spline":
            syst_structs[-1].knobval = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0], dtype=np.float32)
            syst_structs[-1].knob_index = np.array([-3.0, -2.0, -1.0, 0.0, 1.0, 2.0, 3.0], dtype=np.float32)

    for s in syst_structs:
        s.SanityCheck()

    for s in syst_structs:
        s.CreateSpecs(c.m_num_truebins_total if s.mode == "spline" else c.m_num_bins_total)

    prop = profit.PROpeller()
    prop.hist = np.zeros((c.m_num_truebins_total, c.m_num_bins_total))
    prop.histLE = np.zeros((c.m_num_truebins_total,))

    edges = np.array(c.m_channel_truebin_edges)
    edges_lo = edges[:, :-1].flatten()
    edges_hi = edges[:, 1:].flatten()
    prop.histLE[:edges_lo.size] += (edges_lo + edges_hi) / 2.

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("-x", "--xml", type=str, help="XML file name.", default="")
    return parser.parse_args().__dict__

if __name__ == "__main__":
    args = parse_args()
    main(**args)
