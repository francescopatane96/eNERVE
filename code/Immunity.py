
""" Autoimmunity/allergenicity, mouse and conservation modules """

import logging, os, subprocess
from Bio.Blast.Applications import NcbiblastpCommandline
from Bio.Blast import NCBIXML
import pandas as pd
from Protein import *
import shutil

def autoimmunity(list_of_proteins, proteome1, working_dir, NERVE_dir, e_value, minlength, mismatch, substitution) -> list:
        """ Research of human immunogenic and allergenic peptides """

        logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                            filemode = 'a',
                            level = logging.DEBUG,
                            force = True)

        blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir, "database/sapiens_database/sapiens"),
                                            evalue=e_value, outfmt=5, out=os.path.join(working_dir, "sapiens.xml"))
        stdout, stderr = blastx_cline()

        #blast_result = StringIO(stdout)
        for record in NCBIXML.parse(open(os.path.join(working_dir, "sapiens.xml"))):
            query_name = record.query.split(' ')[0]
            for protein in list_of_proteins:
                if query_name in protein.id:
                    tmp_protein = protein
            for alignment in record.alignments:
                for hsp in alignment.hsps:
                    tmp_protein.list_of_shared_human_peps += Protein.hsp_match_parser(hsp.match,
                                                                                    hsp.query,
                                                                                    parsing_window_size=minlength,
                                                                                    max_sub=substitution,
                                                                                    max_mismatch=mismatch)
        os.remove(os.path.join(working_dir, "sapiens.xml"))

        logging.debug('human autoimmunity and allergenicity acquired')
        logging.debug('sum of human peptides starts')

        for protein in list_of_proteins:
            score = 0
            if len(protein.list_of_shared_human_peps) > 0:
                prev_match = protein.list_of_shared_human_peps[0]['match']
                score = len(prev_match)
                for pept in protein.list_of_shared_human_peps[1:]:
                    tmp_match = pept['match']
                    if tmp_match[:len(tmp_match)-1] == prev_match[1:]:
                        score += 1
                    else:
                        score += len(tmp_match)
                    prev_match = tmp_match
            protein.sapiens_peptides_sum = score / protein.length

        mhcpep = pd.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_sapiens.csv"), skipinitialspace=True)
        for protein in list_of_proteins:
            for seq in protein.list_of_shared_human_peps:
                for pep in mhcpep['Epitope.2']:
                    tmp_matchs = Protein.peptide_comparison(seq, pep)
                    protein.list_of_peptides_from_comparison_with_mhcpep_sapiens += tmp_matchs
        return list_of_proteins

def mouse(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, minlength, substitution, mismatch) -> list:
    """ Mouse immunity and allergenicity check """
    logging.basicConfig(filename = os.path.join(working_dir, 'logfile.log'),
                        filemode = 'a',
                        level = logging.DEBUG,
                        force = True)

    blastx_cline = NcbiblastpCommandline(query=proteome1, db=os.path.join(NERVE_dir, "database/mouse_database/mouse"),
                                         evalue=e_value, outfmt=5, out=os.path.join(working_dir, "mouse.xml"))
    stdout, stderr = blastx_cline()
    # outfile = open(os.path.join(working_dir, 'mouse_immunity_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir, "mouse.xml"))):
        query_name = record.query.split(' ')[0]
        tmp_protein = list_of_proteins[0]
        # take the right protein
        for protein in list_of_proteins:
            if query_name in protein.id:  # do not use query_name == p.id
                tmp_protein = protein
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                tmp_protein.list_of_shared_mouse_peps += Protein.hsp_match_parser(hsp.match, hsp.query,
                                                                                          parsing_window_size=minlength,
                                                                                          max_sub=substitution,
                                                                                          max_mismatch=mismatch)
        # print out the peptides (if there are any)
        # if len(tmp_protein.list_of_shared_human_peps) == 0:
        #    outfile.write("\nNo immunogenic peptides for " + query_name)
        # else:
        #    outfile.write("\nList of immunogenic peptides for " + query_name + ": " + str([el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    # outfile.close()
    os.remove(os.path.join(working_dir, "mouse.xml"))  # delete after the computation

    # store peptides from comparison with mouse recognized bacterial mhcpep
    mhcpep = pd.read_csv(os.path.join(NERVE_dir, "database/mhcpep/mhcpep_mouse.csv"), skipinitialspace=True)
    # number_of_proteins = len(list_of_proteins)
    for protein in list_of_proteins:
        for seq in protein.list_of_shared_mouse_peps:
            for pep in mhcpep['Epitope.2']:
                tmp_matches = Protein.peptide_comparison(seq, pep)
                protein.list_of_peptides_from_comparison_with_mhcpep_mouse += tmp_matches
    
    # sum peptides
    logging.debug('mouse sum of peptides starts')
    for protein in list_of_proteins:
        score = 0
        if len(protein.list_of_shared_mouse_peps) > 0:
            prev_match = protein.list_of_shared_mouse_peps[0]['match']
            score = len(prev_match)
            for pept in protein.list_of_shared_mouse_peps[1:]:
                tmp_match = pept['match']
                if tmp_match[:len(tmp_match) - 1] == prev_match[1:]:
                    score += 1
                else:
                    score += len(tmp_match)
                prev_match = tmp_match
        protein.mouse_peptides_sum = score / protein.length
    return list_of_proteins

def conservation(list_of_proteins, working_dir, NERVE_dir, e_value, proteome1, proteome2, minlength,
                 substitution, mismatch) -> list:
    """Conservation analysis between proteome1 and proteome2"""
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)

    bashCmd = f"makeblastdb -in {proteome2} -dbtype prot -parse_seqids -out {os.path.join(working_dir, 'compare_proteome/compare_proteome')}"
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()

    blastx_cline = NcbiblastpCommandline(query=proteome1,
                                         db=os.path.join(working_dir, 'compare_proteome/compare_proteome'),
                                         evalue=e_value, outfmt=5,
                                         out=os.path.join(working_dir, "comparison.xml"))  # 5 is for xml
    stdout, stderr = blastx_cline()

    outfile = open(os.path.join(working_dir, 'conservation_raw_output.txt'), 'w')
    for record in NCBIXML.parse(open(os.path.join(working_dir, "comparison.xml"))):
        query_name = record.query.split(' ')[0]

        for protein in list_of_proteins:
            if query_name in protein.id:  # do not use p.id == query_name
                tmp_protein = protein
                # max_score = 0
        for alignment in record.alignments:
            for hsp in alignment.hsps:
                # if hsp.score > max_score: max_score = hsp.score
                tmp_protein.list_of_shared_conserv_proteome_peps += Protein.hsp_match_parser(hsp.match,
                                                                                                     hsp.query,
                                                                                                     parsing_window_size=minlength,
                                                                                                     max_sub=substitution,
                                                                                                     max_mismatch=mismatch)
        if len(tmp_protein.list_of_shared_human_peps) == 0:
            outfile.write("\nNo shared peptides for " + query_name)
        else:
            outfile.write("\nList of shared peptides for " + query_name + ": " + str(
                [el['match'] for el in tmp_protein.list_of_shared_human_peps]))
    outfile.close()
    # sum peptides
    logging.debug('Run Conservation sum of peptides')
    for protein in list_of_proteins:
        # p.conservation_score = 0
        score = 0
        if len(protein.list_of_shared_conserv_proteome_peps) > 0:
            prev_match = protein.list_of_shared_conserv_proteome_peps[0]['match']
            score = len(prev_match)
            for pept in protein.list_of_shared_conserv_proteome_peps[1:]:
                tmp_match = pept['match']
                if tmp_match[:len(tmp_match) - 1] == prev_match[1:]:
                    score += 1
                else:
                    score += len(tmp_match)
                prev_match = tmp_match
        protein.conservation_score = score / protein.length

    os.remove(os.path.join(working_dir, "comparison.xml"))  # delete after the computation
    if os.path.isdir(os.path.join(working_dir, "compare_proteome")):
        shutil.rmtree(os.path.join(working_dir, "compare_proteome"))
    return list_of_proteins
