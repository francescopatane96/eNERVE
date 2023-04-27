
""" Runs loop razor """

import logging
import os


def razor(list_of_proteins, working_dir, transmem_doms_limit, min_loop_length) -> list:
    " Runs razor module "
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)


    for protein in list_of_proteins:
        if protein.transmembrane_doms >= transmem_doms_limit:
            longest_loop = max(protein.provide_raw_loops(), key=lambda k: len(k), default=None)
            if longest_loop and len(longest_loop) > min_loop_length:
                logging.debug(f'Substituting {str(protein.id)} sequence with its longest loop')
                protein.original_sequence_if_razor = protein.sequence
                protein.sequence = longest_loop
                protein.razored = True
            else:
                logging.debug(f"No replacement found for {str(protein.id)}")
                protein.razored = False
        else:
            protein.razored = False
    return list_of_proteins
'''
    # logging.debug("Warning: razor utilizes X as an exclusive symbol to split the final protein. Check if X is used inside the protein sequences")
    for p in list_of_proteins:
        if p.transmembrane_doms >= transmem_doms_limit:
            longest_loop = max(p.provide_raw_loops(), key = lambda k: len(k))
            if len(longest_loop) > razlen:
                logging.debug(f'Substituting {str(p.id)} sequence with its longest loop')
                p.original_sequence_if_razor = p.sequence
                p.sequence = longest_loop
    return list_of_proteins
'''
