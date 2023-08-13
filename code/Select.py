
""" Protein Selection module """

import pandas as pd
from Protein import *

def select(list_of_proteins, autoimmunity, transmem_doms_limit, padlimit, mouse, antigenlimit, antigen, annotation, razor) -> list:
    """ Selection of suitable vaccine candidate proteins """
    # protein annotations to exclude
    annotations2exclude = ['structural constituent of ribosome', 'DNA binding',
                             'DNA-binding transcription factor activity',
                             'transcription regulator activity', 'rRNA binding', 'RNA binding',
                             'aminoacyl-tRNA ligase activity', 'sequence-specific DNA binding',
                             'catalytic activity, acting on a tRNA', 'catalytic activity, acting on RNA',
                             'tyrosine-tRNA ligase activity', 'aminoacyl-tRNA editing activity',
                             'translation factor activity, RNA binding', 'translation regulator activity',
                             'translation regulator activity, nucleic acid binding',
                             'translation elongation factor activity', 'catalytic activity, acting on DNA'
                             ]
    final_list = []
    
    for p in list_of_proteins:
        # exclude internal proteins with low P_AD or P_ANTIGEN
        if p.localization == "in": continue

        if antigen == "True":
            if p.localization == 'in' and p.p_antigen < antigenlimit and p.p_ad < padlimit: continue
        if antigen != "True":
            if p.localization == 'in' and p.p_ad < padlimit: continue
        if razor == 'False':
            if (p.transmembrane_doms >= transmem_doms_limit) and (p.original_sequence_if_razor is None): continue
        
        if autoimmunity == "True":
            if p.sapiens_peptides_sum > .15: continue
            if len(p.list_of_peptides_from_comparison_with_mhcpep_sapiens) >= 1: continue
            if mouse == "True":
                if p.mouse_peptides_sum > .15: continue
                if len(p.list_of_peptides_from_comparison_with_mhcpep_mouse) >= 1: continue

        annotation_flag = "False"
        if annotation == "True":
            for annot in annotations2exclude:
                if annot in str(p.annotations):
                    annotation_flag = "True"
        if annotation_flag == "True": continue

        final_list.append(p)
    return final_list

def scorer(protein: Protein, mouse: str, autoimmunity: str, antigen: str) -> float:
    """Provides a score for candidate proteins"""

    if mouse == "True" and autoimmunity == "True" and antigen == "True":
        score = (protein.p_ad + protein.p_antigen + \
                    ((protein.reliability_out)) + \
                    (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) + \
                    (1 - (protein.sapiens_peptides_sum / .15)) + \
                    (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)) + \
                    (1 - (protein.mouse_peptides_sum / .15))) / 7
    if mouse != "True" and autoimmunity == "True" and antigen == "True":
        score = (protein.p_ad + protein.p_antigen + \
                    ((protein.reliability_out)) + \
                    (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) + \
                    (1 - (protein.sapiens_peptides_sum / .15))) / 5
    if mouse != "True" and autoimmunity != "True" and antigen == "True":
        score = (protein.p_ad + protein.p_antigen + protein.reliability_out) / 3
    if mouse == "True" and autoimmunity == "True" and antigen != "True":
        score = (protein.p_ad + \
                 ((protein.reliability_out)) + \
                 (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)) + \
                 (1 - (protein.sapiens_peptides_sum / .15)) + \
                 (1 - len(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)) + \
                 (1 - (protein.mouse_peptides_sum / .15))) / 6
    if mouse != "True" and autoimmunity != "True" and antigen != "True":
        score = (protein.p_ad + protein.reliability_out) / 2

    return score

def output(list_of_proteins: list, outfile, mouse: str, autoimmunity: str, antigen: str):
    """Produces output .csv table"""
    df = pd.DataFrame([[str(protein.id),
                        str("".join([str(protein.accession) if protein.accession != None else ""])),
                        (round(scorer(protein, mouse, autoimmunity, antigen), 4)),
                        str(protein.length),
                        str(protein.transmembrane_doms),
                        str(protein.localization),
                        str(protein.reliability_out),
                        # str(", ".join([str(element) for element in protein.localization])),
                        str("".join([str(round(protein.p_antigen, 4)) if protein.p_antigen != None else ""])),
                        str("".join([str(round(protein.p_ad, 4)) if protein.p_ad != None else ""])),
                        str("".join([str(round(protein.conservation_score,
                                                4)) if protein.conservation_score != None else ""])),
                        str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_human_peps if
                                                len(protein.list_of_shared_human_peps) > 0])))),
                        str("".join(str(len([str(dic['match']) for dic in protein.list_of_shared_mouse_peps if
                                                len(protein.list_of_shared_mouse_peps) > 0])))),
                        str("".join(str(len(
                            [str(dic['match']) for dic in protein.list_of_shared_conserv_proteome_peps if
                                len(protein.list_of_shared_conserv_proteome_peps) > 0])))),
                        str("".join([str(round(protein.sapiens_peptides_sum,
                                                   4)) if protein.sapiens_peptides_sum != None else "0"])),
                        str("".join([str(round(protein.mouse_peptides_sum,
                                                   4)) if protein.mouse_peptides_sum != None else "0"])),
                        str("".join([str(protein.annotations) if protein.annotations != None else ""])),
                        str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_sapiens)))),
                        str(", ".join(list(set(protein.list_of_peptides_from_comparison_with_mhcpep_mouse)))),
                        str(protein.sequence),
                        str("".join([str(protein.original_sequence_if_razor) if protein.original_sequence_if_razor != None else ""])),
                        str("".join([str(protein.sequence_out) if protein.sequence_out != None else ""])),
                        str("".join([str(protein.tmhmm_seq) if "M" in str(protein.tmhmm_seq) else ""])),
                        str("".join([str(protein.MHC1_binders) if str(protein.MHC1_binders) != None else ''])),
                        str("".join([str(protein.MHC2_binders) if str(protein.MHC2_binders) != None else ''])),
                        str("".join([str(protein.MHC1_pb_binders) if str(protein.MHC1_pb_binders) != None else ''])),
                        str("".join([str(protein.MHC2_pb_binders) if str(protein.MHC2_pb_binders) != None else ''])),
                        str("".join([str(protein.instability_index) if str(protein.instability_index) != None else ''])),
                        str("".join([str(protein.charge_at_pH_7) if str(protein.charge_at_pH_7) != None else '']))
                       
                        ] for protein in list_of_proteins
                        ],
                        columns=['id',
                                'accession',
                                'score',
                                'length',
                                'transmembrane_doms',
                                'localization',
                                'localization_score',
                                'antigen_probability',
                                'adhesin_probability',
                                'conservation_score',
                                'shared_human_peps',
                                'shared_mouse_peps',
                                'shared_conserved_proteome_peps',
                                'human_peptides_sum',
                                'mouse_peptides_sum',
                                'annotations',
                                'list_of_peptides_from_comparison_with_mhcpep_sapiens',
                                'list_of_peptides_from_comparison_with_mhcpep_mouse',
                                'sequence',
                                'original_sequence_if_razor',
                                'sequence_out',
                                'tmhmm_seq',
                                'MHC1_binders',
                                'MHC2_binders',
                                'MHC1_pb_binders',
                                'MHC2_pb_binders',
                                'Instability_index',
                                'Charge_at_pH7'
                                ]
                        )
    df = df.sort_values(by='score', ascending=False)
    df.to_csv(outfile, index=False)
    
