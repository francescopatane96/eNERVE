
"""Run epitope prediction of proteins"""


from epitopepredict import base, sequtils, analysis, plotting
import math
from utils import *
from Select import *


def epitope(final_proteins, autoimmunity, mouse, mouse_peptides_sum_limit, antigen, working_dir,
            mhci_length, mhcii_length, mhci_overlap, mhcii_overlap, epitope_percentile) -> list:
    """Module to run epitopes prediction"""
    
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    # create predictor object for mhcii
    mhcii_predictor = base.get_predictor('tepitope')
    # create predictor object for mhci
    mhci_predictor = base.get_predictor('basicmhc1')

    # set alleles to consider. Based on literature (Alessandro Sette work) we choose to use supertypes alleles
    m2alleles = base.get_preset_alleles('mhc2_supertypes')
    m1alleles = base.get_preset_alleles('mhc1_supertypes')

    # calculate score for every protein in the list
    protein_scores = []
    #data_i = []
    #data_ii = []
    for protein in final_proteins:
        score = scorer(protein, mouse_peptides_sum_limit, mouse, autoimmunity, antigen)
        protein_scores.append(score)

    if len(protein_scores) != 0:

        # sort scores in ascending order and find selected percentile
        sorted_scores = sorted(protein_scores)
        n = len(sorted_scores)
        percentile_index = math.ceil(n * epitope_percentile) - 1
        percentile = sorted_scores[percentile_index]
	
	
        for p, score in zip(final_proteins, protein_scores):
            if score >= epitope_percentile:
                # run predictions for MHC I and II epitopes
                # use 'threads=0' to use all available cores
                mhci_epitopes = mhci_predictor.predict_sequences(p.sequence, alleles=m1alleles, length=mhci_length,
                                                                 verbose=False, overlap=mhci_overlap)
                mhcii_epitopes = mhcii_predictor.predict_sequences(p.sequence, alleles=m2alleles, length=mhcii_length,
                                                                   verbose=False, overlap=mhcii_overlap)
                
                mhci_epitopes['protein id'] = p.id
                #data_i.append(mhci_epitopes)
                mhcii_epitopes['protein id'] = p.id
                #data_ii.append(mhcii_epitopes)
                
                mhci_epitopes.to_csv('mhci_epitopes_{}.csv'.format(p.id), index=False)
                mhcii_epitopes.to_csv('mhcii_epitopes_{}.csv'.format(p.id), index=False)
                
                p.mhci_epitopes = 'Epitopes acquired. See mhci_epitopes.csv'
                p.mhcii_epitopes = 'Epitopes acquired. See mhcii_epitopes.csv'
                
                #promiscuous binders
                
              
                
    #if len(data_i) > 0:
        #combined_data_i = pd.concat(data_i)
        #combined_data_i = combined_data_i.loc[combined_data_i['rank'] == 1]
        #combined_data_i.to_csv('mhci_epitopes.csv', index=False)
        
    #if len(data_ii) > 0:
        #combined_data_ii = pd.concat(data_ii)
        #combined_data_ii = combined_data_ii.loc[combined_data_ii['rank'] == 1]
        #combined_data_ii.to_csv('mhcii_epitopes.csv', index=False)
                
    return final_proteins
