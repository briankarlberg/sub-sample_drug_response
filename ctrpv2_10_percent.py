# Deploy CTRPv2 sample response on Exacloud
# from arg_test.py on local
# call with: bash hard_path_local.sh <drg_ndx_int>
# then on exa: sbatch hard_path.sh
    # with drug index range in SBATCH array
    
# put data/coding_genes_CCLE.tsv for final exa script - done

print('2022-07-20 standard output check')
import argparse

parser = argparse.ArgumentParser(description = 'CTRPv2 sample response args')
parser.add_argument('drug_index', type = int, help='Input index postion of drug')

args = parser.parse_args() # * Turn off for devel

drug_index = args.drug_index # from local cmnd line / sbatch array on exa
# drug_index = 3 # notebook interactive mode

date = '2022-07-21' # script created

import pandas as pd
from sklearn.feature_selection import RFE
from sklearn import tree
from sklearn.tree import DecisionTreeRegressor
import time

CTRPv2_dups_aved = pd.read_csv('data/CTRPv2_dups_aved.tsv',
                       sep = '\t',
                              index_col = 0)

CTRPv2_dups_aved_uniq_drgs = CTRPv2_dups_aved.Drug.unique()

# Read-in the full coding genes, on for production, 15 seconds i/o
# print('starting expression file read')
read_time_strt = time.time()
raw_feats = pd.read_csv(
    # 'data/coding_genes_CCLE.tsv', # exa
    'data/coding_genes_CCLE.tsv',
    sep = '\t', index_col = 0
)
# print('read time')
# print(time.time() - read_time_strt)
# print(' ')

clf = tree.DecisionTreeRegressor() # classifier for RFE

selector = RFE(clf, n_features_to_select=10, step=.5) # Toogle step percent cut

drg = CTRPv2_dups_aved_uniq_drgs[drug_index]

# Frame of all cell line responses for the this drug    
drug_sub = CTRPv2_dups_aved[CTRPv2_dups_aved.Drug == drg]
# print(drg)
# print('cell line count ' + str(len(drug_sub)))

# subset on only X_expression values for this drug
# then delete full CCLE object to save memory
# X_exp = coding_genes_CCLE[coding_genes_CCLE.index.isin(drug_sub.Cell_line)] # production on

# pilot_steps = [100]
# pilot_steps = [100, 200]
# for i in tqdm(pilot_steps): # Inspection
# for smp_sz in tqdm([750]): # Sample step sizes first, resample loop next

# selector_out = [] # List to direct capture raw RFE output
step_time_strt = time.time()

# for smp_sz in [250]:
for smp_sz in list(range(50,len(drug_sub),50)): # sample step sizes, auto len, for production  
    sel_feat_dct = {}
    for rs in list(range(0,30)):
    # for rs in list(range(0,1)): # Devel
    # for rs in tqdm(list(range(0,30))): # Inspection
        # print('Resampling ' + str(rs))
        rs_name = 'rs_'+str(rs).zfill(2)
        rs_loop_strt = time.time()
        sub_samp_drug_sub = drug_sub.sample(n=smp_sz)
        X_exp = raw_feats[raw_feats.index.isin(sub_samp_drug_sub.Cell_line)]
        # print(len(sub_samp_drug_sub))

        y = sub_samp_drug_sub.Response
        # start_time = time.time()

        selector = selector.fit(X_exp, y) # Fit RFE to the full expression and response
        sel_feat_lst = list(raw_feats.columns[selector.get_support(1)]) # For bio plots

        # res_tup = (drg, selector.support_, selector.ranking_) # Direct caputure raw RFE out
        # selector_out.append(res_tup)

        featuresDF = raw_feats.loc[:, selector.support_] # Keep on, clf built right in now

#         featuresDF.to_csv('results/featuresDF_'+str( # Turn off, re-subset locally
                                                                    # if needed for add clf
#             drg)+'_'+str(
#             smp_sz)+'_'+str(
#                 round(time.time() - loop_strt,1))+'_.tsv', sep = '\t')
        
        smp_len = len(sub_samp_drug_sub)

        X_len = round(smp_len * .8)
        
        trn_tst_splt = {}
        for splt in list(range(0,10)): # Train test split
        # for splt in tqdm(list(range(0,10))): # Train test split
            X = sub_samp_drug_sub.sample(X_len) # training samples
            
            X_trn = featuresDF[featuresDF.index.isin(X.Cell_line)]
            y = sub_samp_drug_sub[sub_samp_drug_sub.index.isin(X.index)].Response
            
            regr_1 = DecisionTreeRegressor(max_depth=2)
            # regr_2 = DecisionTreeRegressor(max_depth=5)
            regr_1.fit(X_trn,y)

            X_tst = featuresDF.loc[
                sub_samp_drug_sub[~sub_samp_drug_sub.Cell_line.isin(X.Cell_line)].Cell_line,
                    :]

            predictions = regr_1.predict(X_tst)
            splt_name = 'splt_'+str(splt)
            
            # trn_tst_splt[splt_name] = predictions
            trn_tst_splt[splt_name] = list(zip(list(X_tst.index),
                                                     predictions))
            # break # train test split

        predict_out = pd.DataFrame(trn_tst_splt) # For accuracy plots
    
        predict_out.to_csv('results/CTRPv2_ten_percent_exa/prd_'+str(
                drg)+'_'+str(
                smp_sz).zfill(3)+'_'+str(rs).zfill(2)+'_'+str( # smp_sz is the sample size,
                                                        # rs is resampling index
                    round(time.time() - rs_loop_strt)).zfill(4)+'_.tsv',
                           sep = '\t')
        
        sel_feat_dct[rs_name] = sel_feat_lst
        # break # rs loop
    feats_out = pd.DataFrame(sel_feat_dct)
        
    feats_out.to_csv('results/CTRPv2_ten_percent_exa/fts_'+str(
                drg)+'_'+str(
                smp_sz).zfill(3)+'_'+str(rs).zfill(2)+'_'+str( # sm is the sample size, rs is resample
                    round(time.time() - step_time_strt)).zfill(4)+'_.tsv',
                           sep = '\t')
