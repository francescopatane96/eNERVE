
import argparse
import numpy as np
import os
from Bio import SeqIO
import logging
from Bio.Seq import Seq
from Protein import *
import subprocess

def load_fasta_file(fasta_file):
    sequences = []
    sequence_ids = []
    with open(fasta_file) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                sequence_ids.append(line[1:])
            else:
                sequences.append(line)
    return sequence_ids, sequences


def extract_features(list_of_proteins, NERVE_dir, iFeature_dir, working_dir, proteome1) -> list:
    """Extract features from protein sequences with iFeature and store them as numpy array in Protein.model_raw_data for subsequent preprocessing and model prediction. These features will be used by adhesin and virulent factor predictos.
    param: working_dir: str
    param: iFeauture_dir: str
    param: proteome1: str
    output: list_of_proteins: pupulates model_raw_data in every Protein element and returns updated list_of_proteins
    """
    #
    features = ["AAC", "DPC", "DDE", "GAAC", "GDPC", "GTPC", "CTDC", "CTDT", "CTDD"]
    features1 = ["AAC", "DPC", "CTDC", "CTDT", "CTDD"]
    extension = ".out"

    # run iFeature
    for feature in features:
        bashCmdMethod(f"python3 {iFeature_dir}/iFeature.py --file {proteome1} --type {feature}\
        --out {os.path.join(working_dir, feature + '.out')}")
        # parse files and update Protein entires
        datasets = [[] for feature in features1]
    for i in range(len(features1)):
        with open(os.path.join(working_dir, features1[i] + extension)) as f:
            lines = f.readlines()[1:]
            for line in lines:
                information = line.split('\t')
                for protein in list_of_proteins:
                    if information[0] in protein.id:
                        protein.model_raw_data1.append(np.array([float(el) for el in information[1:]]))
    # parse files and update Protein entires
    datasets = [[] for feature in features]
    for i in range(len(features)):
        with open(os.path.join(working_dir, features[i] + extension)) as f:
            lines = f.readlines()[1:]
            for line in lines:
                information = line.split('\t')
                for protein in list_of_proteins:
                    if information[0] in protein.id:
                        protein.model_raw_data.append(np.array([float(el) for el in information[1:]]))
    # delete files after computation
    for file in features:
        os.remove(os.path.join(working_dir, file + extension))
    # standardization
    for protein in list_of_proteins:
        # concatenate
        protein.model_raw_data = np.concatenate(protein.model_raw_data)
    for protein in list_of_proteins:
        # concatenate
        protein.model_raw_data1 = np.concatenate(protein.model_raw_data1)

    return list_of_proteins

def bashCmdMethod(bashCmd):
    """Run bash commands
    param: bashCmd: bash command to be run"""
    process = subprocess.Popen(bashCmd.split(), stdout=subprocess.PIPE)
    output, error = process.communicate()
    return output, error

class protein_element:
    """Class to handle fasta file elements similarly to the biopython fasta file parser"""
    def __init__(self, name, seq):
        self.name=name
        self.seq=seq

def proteome_downloader(working_dir, proteome_id, filename='input_proteome.fasta', output_dir=os.getcwd(), format_ = "fasta") -> None:
    """Downloads proteome from uniprot database into output multi-fasta file
    param: proteome_id: uniprot unique proteome id, not case-sensitive
    param: output_dir: output directory (default: current directory)
    param: format: uniprot API required format (default:fasta)
    param: filename: output proteome filename (default: input_proteome.fasta)
    """
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)
    try:
        url=f'https://rest.uniprot.org/uniprotkb/stream?compressed=false&format={format_}&query=%28proteome%3A{proteome_id}%29'
        response = requests.get(url, stream = True)
        text_file = open(os.path.join(output_dir, filename), 'wb')
        for chunk in response.iter_content(chunk_size=1024):
              text_file.write(chunk)
        # raise an AssertionError if the given proteome ID is not valid
        assert text_file.tell() > 0, 'First download attempt failed'
        text_file.close()
    except AssertionError:
        logging.debug(f'28proteome is not a valid building block')
        try:
            url=f'https://rest.uniprot.org/uniprotkb/stream?format={format_}&query=%28proteome%3A{proteome_id}%29'
            response = requests.get(url, stream = True)
            text_file = open(os.path.join(output_dir, filename), 'wb')
            for chunk in response.iter_content(chunk_size=1024):
                  text_file.write(chunk)
            # raise an AssertionError if the given proteome ID is not valid
            assert text_file.tell() > 0, 'Second download attempt failed'
            text_file.close()
        except AssertionError:
            logging.debug(f'avoiding compression is not a valid building block')
            try:
                url=f'https://rest.uniprot.org/uniparc/stream?format={format_}&query=%28upid%3A{proteome_id}%29'
                response = requests.get(url, stream = True)
                text_file = open(os.path.join(output_dir, filename), 'wb')
                for chunk in response.iter_content(chunk_size=1024):
                      text_file.write(chunk)
                # raise an AssertionError if the given proteome ID is not valid
                assert text_file.tell() > 0, 'Third download attempt failed'
                text_file.close()
            except AssertionError as e:
                logging.debug(f'28upid is not a valid building block')
                print(f'Unable to download proteome {proteome_id} due to invalid proteome ID, Uniprot API failure or wrong path. In case of uniprot API failure provide the proteome as file.')
                raise SystemExit(e)
    return None

def proteome_uploader(infile:str)->list:
    """Function to read and parse fasta files. Bio SeqIO is not suitable because it chops sequence names. It
    will be used only to validate fasta file format.
    param: infile: path to fasta file"""
    class protein_element:
        """Class to handle fasta file elements similarly to the biopython fasta file parser"""
        def __init__(self, name, seq):
            self.name=name
            self.seq=seq

    proteome_elements = []
    proteome_data = {}
    name = ''  # initialize name variable
    infile = open(infile, 'r').readlines()
    for i in range(len(infile)):
        if infile[i].startswith('>'):
            name = infile[i].strip()[1:]
            proteome_data[name] = ''
        elif name != '':
            proteome_data[name] += infile[i].strip()
    for element in proteome_data:
        proteome_elements.append(protein_element(element, proteome_data[element]))
    return proteome_elements


def is_fasta(filename:str):
    """Function that rise an error if the format is not .fasta.
    param: filename: path to fasta file"""
    with open(filename, "r") as handle:
        fasta = list(SeqIO.parse(handle, "fasta"))
        # biopython silently fails if the format is not fasta returning an empty generator
        # any() returns False if the list is empty
        if any(fasta) == True:
            fasta=proteome_uploader(filename)
            return fasta
        else:
            raise ValueError(f'{filename} is not in fasta format')

def dir_path(path:str)->str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path

def quality_control(path_to_fasta: str, working_dir: str, upload=False) -> dir_path:
    """
    Remove sequences with non-canonical aminoacid symbols. U (Se-Cys) is substituted with C (Cys). Returns
    {working_dir}/discarded_sequences_{input_fasta_file_name}.fasta containing discarded sequences and\
    {working_dir}/cleaned_{input_fasta_file_name}.fasta with cleaned sequences
    param: path_to_fasta: full path to fasta file containing the proteome with .fasta extension. Input fasta file\
    will not be overwritten
    param: working_dir: working directory, were cleaned fasta file will be saved
    param: upload: if True the function returns a list containing the filtered sequences as protein_element objects
    output: path to the cleaned fasta file named as {working_dir}/cleaned_{input_fasta_file_name}.fasta
    """
    # define logging file
    logging.basicConfig(filename=os.path.join(working_dir, 'logfile.log'),
                        filemode='a',
                        level=logging.DEBUG,
                        force=True)

    aa_dic = {'C': 'C', 'D': 'D', 'S': 'S', 'Q': 'Q', 'K': 'K', 'I': 'I', 'P': 'P', 'T': 'T', 'F': 'F', 'N': 'N',
              'G': 'G', 'H': 'H', 'L': 'L', 'R': 'R', 'W': 'W', 'A': 'A', 'V': 'V', 'E': 'E', 'Y': 'Y', 'M': 'M',
              'U': 'C'}
    filtered_sequences, discarded_sequences = [], []
    # control formatting
    fasta_list = is_fasta(path_to_fasta)
    # filename needed to create the output file
    file_name = os.path.basename(path_to_fasta)
    output_file = os.path.join(working_dir, "_".join(["cleaned", file_name]))
    output_discarded_sequences = os.path.join(working_dir, "_".join(["discarded_sequences", file_name]))
    for record in fasta_list:
        flag = True
        new_seq = ''
        # check sequence
        for aa in str(record.seq):
            if aa not in aa_dic:
                flag = False
                logging.debug(f'Found non-canonical aminoacid "{aa}" in sequence: {record.name}')
            elif aa == "U":
                logging.debug(
                    f'Found non-canonical aminoacid "{aa}" (Selenocysteine) in sequence: {record.name}, substituting to Cysteine')
                new_seq += aa_dic[aa]
            else:
                new_seq += aa_dic[aa]
        # check name
        if ">" in record.name:
            logging.debug(f'Found non-canonical character ">" in sequence name:\n{record.name}\nSubstituting with "*"')
            record.name = record.name.replace(">", "*")
        record.seq = Seq(new_seq)

        if flag == True:
            filtered_sequences.append(record)
        else:
            discarded_sequences.append(record)
            logging.debug(f'Sequence {record.name} has been discarded for the presence of non-canonical aminoacids.')
            # output filtered overwriting input fasta file
    filename = open(output_file, 'w')
    for sequence in filtered_sequences:
        filename.write(f'>{str(sequence.name)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    # SeqIO.write(filtered_sequences, filename, "fasta")
    filename.close()
    # output discarded sequences
    filename = open(output_discarded_sequences, 'w')
    for sequence in discarded_sequences:
        filename.write(f'>{str(sequence.name)}\n')
        filename.write(f'{str(sequence.seq)}\n')
    # SeqIO.write(discarded_sequences, filename, "fasta")
    filename.close()
    # repath output_file to be absolute
    output_file = os.path.join('/', os.path.relpath(output_file, start='/'))
    if upload == True:
        # return filtered sequnces and the new path
        return filtered_sequences, output_file
    # return the path to the cleaned file
    return output_file