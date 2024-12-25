import argparse
import profit
import numpy as np
import pandas as pd
import uproot
import awkward as ak

# TODO: WHY DOES SYSTSTRUCT USE unique_ptr's TO CV + SYST UNIVERSES???

def main(xml="", inject=None):
    c = profit.PROconfig(xml)
    print("DONE PARSING")

    # open up the files
    ttree_dfs = []
    for fname, tree in zip(c.m_mcgen_file_name, c.m_mcgen_tree_name):
        ttree_dfs.append(load_df(fname, tree))

    # load systematic weights
    eventweights = []
    for i_f, (f, ttree) in enumerate(zip(c.m_mcgen_file_name, ttree_dfs)):
         evw = loadsysts(f, i_f, ttree, c)
         eventweights.append(evw)

    # map_systematic_num_universe <- number of universes for each systematic
    map_systematic_num_universe = {}
    for col in eventweights[0].columns.get_level_values(0).unique():
         map_systematic_num_universe[col] = len(eventweights[0][col].columns)

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

    print("DONE LOADING SYSTEMATICS")

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
    # TODO: WHY IS THIS THE WRONG SIZE???
    prop.histLE[:edges_lo.size] += (edges_lo + edges_hi) / 2.

    # Process events
    for idf, tdf in enumerate(ttree_dfs):
        process_events(tdf, idf, syst_structs, eventweights[idf], prop, c)

    print("DONE PROCESSING EVENTS")

    # Load systematic weights into PROsyst
    systs = profit.PROsyst(syst_structs)

    osc = profit.PROsc(prop)

    print("DONE LOADING PRO CLASSES")

    # Injected model point
    return

def load_df(fname, treename, branches=None):
    if fname.endswith(".root"):
        tf = uproot.open(fname)
        return make_array(tf[treename], branches)
    elif fname.endswith(".df"):
        df = pd.read_df(fname, key=treename)
        if branches is not None:
            return df[branches]
    else: 
        print(fname)
        assert(False)

def make_array(ttree, branches=None):
    if branches is None:
        branches = ttree.keys()

    return ak.to_dataframe(ttree.arrays(branches, library="ak"), how="outer")

def systematics_df(df, c):
    # Restrict to the configured systematics
    df = df[[k for k in df.columns if k in c.m_mcgen_variation_allowlist and k not in c.m_mcgen_variation_denylist]]
    modes = [c.m_mcgen_variation_type_map[syst] for syst in df.columns]

    def syst_index_name(index, mode):
        if mode == "spline": # Multisigma
            return ["ms3", "ms2", "ms1", "cv", "ps1", "ps2", "ps3"][index]
        else: # Multisim
            return "univ_%i" % index
        # TODO: morph

    # Reformat to the structure we want -- a hierarchical index on the columns, flat index on the rows
    if isinstance(df.index, pd.MultiIndex): 
        assert(df.index.nlevels == 2)
        df = df.unstack()
        df.columns = pd.MultiIndex.from_tuples([(col[0], syst_index_name(col[1], c.m_mcgen_variation_type_map[col[0]])) for col in df.columns])

    # TODO: handle other dataframe formats

    # trim columns with nans
    df = df.dropna(axis=1, how='all')

    return df

def loadsysts(fname, fid, ttree_df, c):
    friends = []
    friend_names = []
    for friend, friend_f in zip(c.m_mcgen_file_friend_treename_map.get(fname, []), c.m_mcgen_file_friend_map.get(fname, [])):
        # Remove duplicates (TODO: why are they there???)
        if friend in friend_names: continue

        friend_names.append(friend)
        friends.append(load_df(friend_f, friend))

    # get the event weights
    # event_weights_pertree = [make_array(t, [k for k in t.keys() if k in c.m_mcgen_variation_allowlist and k not in c.m_mcgen_variation_denylist]) 
    #    for t in [ttree] + friends]
    # event_weights_pertree = [df[[k for k in df.columns if k in c.m_mcgen_variation_allowlist and k not in c.m_mcgen_variation_denylist]]
    #    for df in [ttree_df] + friends]
    event_weights_pertree = [systematics_df(df, c) for df in [ttree_df] + friends]
    event_weights_pertree = pd.concat([df for df in event_weights_pertree if not df.empty], axis=1)

    for ib in range(len(c.m_branch_variables[fid])):
        branch = c.m_branch_variables[fid][ib]
        c.m_branch_variables[fid][ib].branch_formula = profit.DataFrameFormula("branch_form_%i_%i" % (fid, ib), branch.name, ttree_df)
        c.m_branch_variables[fid][ib].branch_true_L_formula = profit.DataFrameFormula("branch_L_form_%i_%i" % (fid, ib), branch.true_L_name, ttree_df)
        c.m_branch_variables[fid][ib].branch_true_value_formula = profit.DataFrameFormula("branch_true_form_%i_%i" % (fid, ib), branch.true_param_name, ttree_df)
        c.m_branch_variables[fid][ib].branch_true_pdg_formula = profit.DataFrameFormula("branch_pdg_form_%i_%i" % (fid, ib), branch.pdg_name, ttree_df)

        if c.m_mcgen_additional_weight_bool[fid][ib]:
            c.m_branch_variables[fid][ib].branch_monte_carlo_weight_formula = profit.DataFrameFormula("branch_add_weight_%i_%i" % (fid, ib), c.m_mcgen_additional_weight_name[fid][ib], ttree_df)

    return event_weights_pertree

def process_branch(c, branch, evws, mcpot, subchannel_index, syst_vector, syst_additional_weight, inprop):
    # Load values
    total_num_sys = len(syst_vector)
    reco_value = branch.GetValue()
    true_param = branch.GetTrueValue()

    baseline = branch.GetTrueL()
    true_value = baseline / true_param
    pdg_id = branch.GetTruePDG()
    run_syst = branch.GetIncludeSystematics()
    mc_weight = branch.GetMonteCarloWeight()
    mc_weight *= c.m_plot_pot / mcpot

    # fix mcweight to be series. TODO: remove fix
    if not isinstance(mc_weight, pd.Series):
        mc_weight = pd.Series(mc_weight, true_param.index)

    model_rule = pd.Series(branch.GetModelRule(), true_param.index)
    run_syst = branch.GetIncludeSystematics()

    global_bin = profit.FindGlobalBin(c, reco_value, subchannel_index)
    global_true_bin = profit.FindGlobalTrueBin(c, true_value, subchannel_index) if run_syst else pd.Series(0, true_param.index)

    valid = (global_bin >= 0) & (global_true_bin >= 0)

    # Fill values in vectors
    inprop.reco = np.concatenate((inprop.reco, reco_value[valid]))
    inprop.added_weights = np.concatenate((inprop.added_weights, mc_weight[valid]))
    inprop.bin_indices = np.concatenate((inprop.bin_indices, global_bin[valid]))
    inprop.pdg = np.concatenate((inprop.pdg, pdg_id[valid]))
    inprop.truth = np.concatenate((inprop.truth, true_param[valid]))
    inprop.baseline = np.concatenate((inprop.baseline, baseline[valid]))
    inprop.model_rule = np.concatenate((inprop.model_rule, model_rule[valid]))
    inprop.true_bin_indices = np.concatenate((inprop.true_bin_indices, global_true_bin[valid]))
    inprop.hist[global_true_bin[valid], global_bin[valid]] += mc_weight[valid]

    if not run_syst: return

    # Systematics
    for i,s in enumerate(syst_vector):
        additional_weight = syst_additional_weight[i]
        evw = evws[s.GetSysName()]
        if s.mode == "spline":
            s.FillCV(global_true_bin[valid], mc_weight[valid])
            for i_univ in range(s.GetNUniverse()):
                col = evw.columns[i_univ]
                s.FillUniverse(i_univ, global_true_bin[valid], (mc_weight*additional_weight*evw[col])[valid])
        # TODO: WHY DO NON-SPLINE SYSTEMATICS USE RECO VARIABLE???
        else:
            s.FillCV(global_bin[valid], mc_weight[valid])
            for i_univ in range(s.GetNUniverse()):
                col = evw.columns[i_univ]
                s.FillUniverse(i_univ, global_bin[valid], (mc_weight*additional_weight*evw[col])[valid])
    # Done!
    return

def process_events(df, fid, syst_structs, evw, prop, c):
    sys_weight_formula = []
    for s in syst_structs:
        if s.HasWeightFormula():
            sys_weight_formula.append(profit.DataFrameFormula("weightMapsFormulas_%i_%s" % (fid, s.GetSysName()), s.GetWeightFormula(), df))

    branches = c.m_branch_variables[fid]
    subchannel_index = []
    for ib, b in enumerate(branches):
        subchannel_name = c.m_branch_variables[fid][ib].associated_hist
        subchannel_index.append(c.GetSubchannelIndex(subchannel_name))

    # Calculate the sys weights
    sys_weight_values = []
    for s in sys_weight_formula:
        sys_weight_values.append(s.EvalInstance())

    # calculate stuff on each event
    for ib, b in enumerate(branches):
        process_branch(c, b, evw, c.m_mcgen_pot[fid], subchannel_index[ib], syst_structs, sys_weight_values, prop)

    # Done!
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
