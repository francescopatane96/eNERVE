
from Bio.SeqUtils.ProtParam import ProteinAnalysis


def protparam(list_of_proteins):
    
    for p in list_of_proteins:
        analysed_seq = ProteinAnalysis(p.original_sequence_if_razor if p.original_sequence_if_razor != None else p.sequence)
        p.instability_index = analysed_seq.instability_index()
        p.charge_at_pH_7 = analysed_seq.charge_at_pH(7)
        
    return list_of_proteins
