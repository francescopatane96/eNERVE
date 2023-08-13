
# Author: Francesco Patanè, University of Padova, SynBio lab (Supervisor: Prof. Francesco Filippini)

"""Run eNERVE 1.0, eukaryotic New Enhanced Reverse Vaccinology Environment"""


def warn(*args, **kwargs):
    pass
from signal import signal, SIGPIPE, SIG_DFL
signal(SIGPIPE, SIG_DFL)
import warnings
warnings.warn = warn
import os
import logging
os.environ['TF_CPP_MIN_LOG_LEVEL'] = '3' 
import tensorflow
logging.basicConfig(level=logging.WARNING)
import time
from Subcellular import *
from Adhesin import *
from Antigenicity import *
from Topology import *
from Razor import *
from Function import *
from Immunity import *
from Epitope import *
from Protein import *
import requests
from Protparam import *


def dir_path(path:str) -> str:
    '''Path validator'''
    if os.path.isdir(path) == False:
        raise argparse.ArgumentTypeError(f'{path} is not a valid path')
    return path

class Args(NamedTuple):
    '''Command-line arguments'''

    annotation:str
    autoimmunity:str
    topology:str
    e_value:float
    minlength:int
    mismatch:int
    mouse:str
    proteome1:str
    proteome2:str
    padlimit:float
    razor:str
    antigenlimit:float
    loclimit:float
    min_loop_length:int
    select:str
    substitution:float
    transmem_doms_limit:int
    antigen:str
    working_dir:str
    NERVE_dir:str
    iFeature_dir:str
    DeepFri_dir:str
    epitopes:str
    mhci_length:int
    mhcii_length:int
    mhci_overlap:int
    mhcii_overlap:int
    epitope_percentile:float


    def print_args(self):
        return (f'''
                    annotation: {self.annotation},
                    autoimmunity: {self.autoimmunity},
                    topology: {self.topology},
                    e_value: {self.e_value},
                    minlength: {self.minlength}, 
                    mismatch: {self.mismatch},
                    mouse: {self.mouse},
                    proteome1: {self.proteome1},
                    proteome2: {self.proteome2},
                    padlimit: {self.padlimit},
                    razor: {self.razor},
                    antigenlimit: {self.antigenlimit},
                    loclimit: {self.loclimit},
                    min_loop_length_razorlen: {self.min_loop_length},
                    select: {self.select},
                    substitution: {self.substitution},
                    transmem_doms_limit: {self.transmem_doms_limit},
                    antigen: {self.antigen},
                    working_dir: {self.working_dir},
                    NERVE_dir: {self.NERVE_dir},
                    iFeature_dir: {self.iFeature_dir},
                    DeepFri_dir: {self.DeepFri_dir},
                    epitopes: {self.epitopes},
                    mhci_length: {self.mhci_length},
                    mhcii_length: {self.mhcii_length},
                    mhci_overlap: {self.mhci_overlap},
                    mhcii_overlap: {self.mhcii_overlap},
                    epitope_percentile: {self.epitope_percentile}
                    ''')

def get_args() -> Args:
    '''Get command-line arguments'''

    parser = argparse.ArgumentParser(
        description = "Run vaccine candidate prediction",
        formatter_class = argparse.ArgumentDefaultsHelpFormatter)

    parser.add_argument('-a', '--annotation',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of annotation module. Uses DeepFri to retrieve protein functional onthologies.",
                        type=str,
                        default="True",
                        required=False
                        )
    parser.add_argument('-ai', '--autoimmunity',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of autoimmunity module.",
                        type=str,
                        default="True",
                        required=False
                        )
    parser.add_argument('-tp', '--topology',
                        metavar='\b',
                        help="Topology predictor to use: TMHMM or deepTMHMM",
                        choices=['tmhmm', 'deeptmhmm'],
                        default='tmhmm',
                        required=False
                        )
    parser.add_argument('-ev', '--e_value',
                        metavar='\b',
                        help="Expect-value used in blastp for auto-immunity module",
                        type=float,
                        default=1e-10,
                        required=False,
                        )
    parser.add_argument('-ml', '--minlength',
                        metavar='\b',
                        help="Minimal length required for shared peptides to be extracted in comparison analyses versus human and/or mouse",
                        type=int,
                        default=9,
                        required=False,
                        )
    parser.add_argument('-mm', '--mismatch',
                        metavar='\b',
                        help="Maximal number of not compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules",
                        type=int,
                        default=1,
                        required=False,
                        )
    parser.add_argument('-m', '--mouse',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of the mouse immunity module. This module compares proteome1 with mouse proteome and a further analysis of the eventual shared peptides is carried out as in the autoimmunity module",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-p1', '--proteome1',
                        metavar='\b',
                        help='Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=True,
                        )
    parser.add_argument('-p2', '--proteome2',
                        metavar='\b',
                        help='Path to proteome or Uniprot proteome ID (see: https://www.uniprot.org/proteomes/?query=&sort=score).',
                        type=str,
                        required=False,
                        )
    parser.add_argument('-pl', '--padlimit',
                        metavar='\b',
                        help="Set the probability of adhesin (pad) value cut-off for proteins with 'in' localization in the select module. Thus, these proteins with a pad value < cut-off are discarded (0.-1)",
                        type=float,
                        default=0.51,
                        required=False,
                        )
    parser.add_argument('-rz', '--razor',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of the loop-razor module. This module allows the recovery of protein vaccine candidates with more than 0 transmembrane domains, that would otherwise be discarded in the select module. The longest loop with minimum len == 'razlen' aa will replace the original protein sequence for following NERVE steps",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-ig', '--antigenlimit',
                        metavar='\b',
                        help="Cut-off value for euANTIGEN in the select module (0.-1)",
                        type=float,
                        default=0.80,
                        required=False,
                        )
    parser.add_argument('-locl', '--loclimit',
                        metavar='\b',
                        help="Set the probability of localization (ploc) value cut-off for proteins with 'in' localization in the select module. Thus, these proteins with a pad value < cut-off are discarded (0.-1)",
                        type=float,
                        default=0.60,
                        required=False,
                        )
    parser.add_argument('-rl', '--min_loop_length',
                        metavar='\b',
                        help="Set minimal length of loop considered in loop-razor module",
                        type=int,
                        default=9,
                        required=False,
                        )
    parser.add_argument('-s', '--select',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of select module, which filters PVC from proteome1",
                        type=str,
                        default="True",
                        required=False,
                        )
    parser.add_argument('-ss', '--substitution',
                        metavar='\b',
                        help="Maximal number of compatible substitutions allowed in shared peptides alignment windows of 'minlength' size in immunity modules",
                        type=int,
                        default=3,
                        required=False,
                        )
    parser.add_argument('-tdl', '--transmem_doms_limit',
                        metavar='\b',
                        help="Parameter of select module. Proteins with trasmembrane domains >= transmem_doms_limit are discarded",
                        type=int,
                        default=2,
                        required=False,
                        )
    parser.add_argument('-ang', '--antigen',
                        metavar='\b',
                        help="Activation (True) or deactivation (False) of euANTIGEN module, predictor of the probability of being an antigen",
                        type=str,
                        default="False",
                        required=False,
                        )
    parser.add_argument('-wd', '--working_dir',
                        metavar='\b',
                        help='path to working directory. If not existing, a working directory with the given path is created',
                        type=str,
                        required=False,
                        default='./'
                        )
    parser.add_argument('-nd', '--NERVE_dir',
                        metavar='\b',
                        help='path to NERVE folder',
                        type=dir_path,
                        required=False,
                        default='./'
                        )
    parser.add_argument('-id', '--iFeature_dir',
                        metavar='\b',
                        help='NERVE folder',
                        type=dir_path,
                        required=False,
                        default='./iFeature'
                        )
    parser.add_argument('-dfd', '--DeepFri_dir',
                        metavar='\b',
                        help='Path to DeepFri folder (download from: https://github.com/flatironinstitute/DeepFRI)',
                        type=dir_path,
                        required=False,
                        default='./DeepFRI'
                        )
    parser.add_argument('-ep', '--epitopes',
                        metavar='\b',
                        type=str,
                        help='Activate or deactivate epitopes module',
                        required=False,
                        default="True"
                        )
    parser.add_argument('-m1l', '--mhci_length',
                        metavar='\b',
                        type=int,
                        help='mhci binders length (9, 10, 11 are available)',
                        required=False,
                        choices=[9, 10, 11],
                        default=9
                        )
    parser.add_argument('-m2l', '--mhcii_length',
                        metavar='\b',
                        type=int,
                        help='mhcii binders length (9, 11, 12, 15 are available)',
                        required=False,
                        choices=[9, 11, 13, 15],
                        default=11
                        )
    parser.add_argument('-m1ovr', '--mhci_overlap',
                        metavar='\b',
                        type=int,
                        help='mhci-epitope overlap',
                        required=False,
                        choices=[1, 2],
                        default=1
                        )
    parser.add_argument('-m2ovr', '--mhcii_overlap',
                        metavar='\b',
                        type=int,
                        help='mhcii-epitope overlap',
                        required=False,
                        choices=[1, 2],
                        default=1
                        )
    parser.add_argument('-prt', '--epitope_percentile',
                        metavar='\b',
                        type=float,
                        help='percentile decision threshold on whick to predict epitopes from full length proteins',
                        required=False,
                        default=0.75
                        )
                        
                          
    args = parser.parse_args()

    return Args(args.annotation, args.autoimmunity, args.topology, args.e_value, args.minlength, args.mismatch,
                args.mouse, args.proteome1, args.proteome2, args.padlimit, args.razor, args.antigenlimit, args.loclimit,
                args.min_loop_length, args.select, args.substitution, args.transmem_doms_limit,
                args.antigen, args.working_dir, args.NERVE_dir, args.iFeature_dir, args.DeepFri_dir, args.epitopes,
                args.mhci_length, args.mhcii_length, args.mhci_overlap, args.mhcii_overlap, args.epitope_percentile)

def main():
    """ Run eNERVE """
    # record time:
    enerve_start = time.time()
    args = get_args()

    header = 'eNERVE 1.0 starts. [Created by Francesco Patanè under the supervision of Prof. Francesco Filippini, University of Padova]'
    print("=" * 100)
    print("\n{:^50}\n".format(header))
    print("=" * 100)

    # init workdir:
    if args.working_dir[-1] != '/':
        args = args._replace(working_dir = args.working_dir + '/')
    # create working dir if it doesn't exist
    if os.path.isdir(args.working_dir)  == False:
        os.makedirs(args.working_dir)

    for handler in logging.root.handlers[:]:
        logging.root.removeHandler(handler)
    # define log file
    logging.basicConfig(filename = os.path.join(args.working_dir, 'logfile.log'),
                        filemode = 'w',
                        level = logging.DEBUG)
    logging.debug(f'Running eNERVE with the following parameters: \n {args.print_args()}')
    print(f'Running eNERVE with the following parameters: \n {args.print_args()}')

    # 1.check input and download proteome
    if os.path.isfile(args.proteome1) == True:
        logging.debug(f'{args.proteome1} found as {args.proteome1}')
        # repath proteome as absolute path
        args = args._replace(proteome1 = os.path.join('/', os.path.relpath(args.proteome1, start='/')))
    elif os.path.isfile(os.path.join(args.working_dir, args.proteome1)) == True:
        logging.debug(f'{args.proteome1} was found in {args.working_dir}')
        args = args._replace(proteome1 = os.path.join(args.working_dir, args.proteome1))
    else:
        logging.debug(f'{args.proteome1} is not a file, download from Uniprot.')
        try:
            proteome_downloader(args.working_dir, args.proteome1, filename=os.path.join(args.working_dir, \
                                                                                        'proteome1.fasta'))
        except Exception as e:
            raise ValueError(f'{args.proteome1} raised the following error:\n{e}')
        logging.debug(f'{args.proteome1} successfully downloaded')
        args = args._replace(proteome1=os.path.join(args.working_dir, 'proteome1.fasta'))

    if args.proteome2:
        if os.path.isfile(args.proteome2) == True:
            logging.debug(f'{args.proteome2} found as {args.proteome2}')
            # repath proteome as absolute path
            args = args._replace(proteome2=os.path.join('/', os.path.relpath(args.proteome2, start='/')))
        elif os.path.isfile(os.path.join(args.working_dir, args.proteome2)) == True:
            logging.debug(f'{args.proteome2} was found in {args.working_dir}')
            args = args._replace(proteome2=os.path.join(args.working_dir, args.proteome2))
        else:
            logging.debug(f'{args.proteome2} is not a file, download from Uniprot.')
            try:
                proteome_downloader(args.working_dir, args.proteome2, filename=os.path.join(args.working_dir, \
                                                                                            'proteome2.fasta'))
            except:
                raise logging.error(f'{args.proteome2} raised the following error:\n{e}')
            logging.debug(f'{args.proteome2} successfully downloaded')
            args = args._replace(proteome2=os.path.join(args.working_dir, 'proteome2.fasta'))

    print("=" * 50)
    print("{:^50}".format('Quality control starts'))
    print("=" * 50)

    # 2.Quality control
    start = time.time()
    logging.debug(f'Start quality control of proteome1 ({args.proteome1})')
    # during the quality control, upload sequences from proteome1
    list_of_fasta_proteins, proteome1_new_path = quality_control(args.proteome1, args.working_dir, upload=True)
    # update input path of proteome1
    args = args._replace(proteome1=proteome1_new_path)
    logging.debug(f'Finish quality control of proteome1. Updated path: ({args.proteome1})')
    if len(list_of_fasta_proteins) == 0:
        raise ValueError(
            f'All input protein sequences have been discarded. See {os.path.join(args.working_dir, "logfile.log")} for more information.')
    if args.proteome2:
        logging.debug(f'Start quality control of proteome2 ({args.proteome2})')
        proteome2_new_path = quality_control(args.proteome2, args.working_dir)
        # update proteome2 new path
        args = args._replace(proteome2=proteome2_new_path)
        logging.debug(f'Finish quality control of proteome2. Updated path: ({args.proteome2})')
    logging.debug(f'Extract protein sequences and IDs from proteome1')
    list_of_proteins = []
    for p in list_of_fasta_proteins:
        p_id = str(p.name)
        p_seq = str(p.seq)
        list_of_proteins.append(Protein(p_id, p_seq))
    end = time.time()
    logging.debug(f'{len(list_of_fasta_proteins)} proteins loaded in {end - start} seconds')

    print(f'{len(list_of_fasta_proteins)} proteins loaded')

    

    print("=" * 50)
    print("{:^50}".format('Descriptors calculation for built-in ML models starts'))
    print("=" * 50)

    start = time.time()
    list_of_proteins = extract_features(list_of_proteins, args.NERVE_dir, args.iFeature_dir, args.working_dir,
                                        args.proteome1)
    end = time.time()

    print("=" * 50)
    print("{:^50}".format('Finish descriptors calculation, predictions start'))
    print("=" * 50)

    logging.debug(f'Finish descriptors calculator in {end - start} seconds. Saved in {args.proteome1}')

    logging.debug(f'Start predictions')

    print("=" * 50)
    print("{:^50}".format('Making predictions, euLOCALIZATION starts'))
    print("=" * 50)

    # 3.Subcellular localization module
    start = time.time()
    logging.debug('sucellular localization starts with ...')
    list_of_proteins = euloc(list_of_proteins, args.working_dir, args.NERVE_dir, args.loclimit)
    end = time.time()
    logging.debug(f'Done run in {end - start} seconds')

    print("=" * 50)
    print("{:^50}".format('Localization acquired, euSPAAN starts'))
    print("=" * 50)

    # 4.Adhesins and adhesins-like module
    start = time.time()
    logging.debug('Adhesins prediction starts with euSPAAN')
    list_of_proteins = euspaan(list_of_proteins, args.working_dir, args.NERVE_dir)
    end = time.time()
    logging.debug(f'euSPAAN done run in {end - start} seconds')

    print("=" * 50)
    print("{:^50}".format('Adhesins acquired'))
    print("=" * 50)

    # 5.Antigenicity module
    if args.antigen == "True":
        print("=" * 50)
        print("{:^50}".format('euANTIGEN starts'))
        print("=" * 50)
        start = time.time()
        logging.debug('Antigenicity prediction with euANTIGEN starts')
        list_of_proteins = euntigen(list_of_proteins, args.working_dir, args.NERVE_dir)
        end = time.time()
        logging.debug(f'Done run in {end - start} seconds')

        print("=" * 50)
        print("{:^50}".format('Antigens acquired'))
        print("=" * 50)

   # 6.TMhelices/topology module
    logging.debug('Topology prediction starts')

    if args.topology == 'tmhmm':
        print("=" * 50)
        print("{:^50}".format('Topology with TMHMM starts'))
        print("=" * 50)
        start = time.time()
        list_of_proteins = tmhelices(list_of_proteins, args.working_dir)
        end = time.time()
        logging.debug(f'TMHMM done in {end - start} seconds')

        print("=" * 50)
        print("{:^50}".format('Topology acquired'))
        print("=" * 50)

    elif args.topology == 'deeptmhmm':
        start = time.time()
        list_of_proteins = tmhelices_deep(list_of_proteins, args.working_dir)
        end = time.time()
        logging.debug(f' deepTMHMM done in {end - start} seconds')
        print('topology acquired')

    # 7.Razor module
    if args.razor == "True":
        print("=" * 50)
        print("{:^50}".format('loops with Razor starts'))
        print("=" * 50)
        logging.debug("Loop-razor starts ...")
        start = time.time()
        list_of_proteins = razor(list_of_proteins, args.working_dir, args.transmem_doms_limit, args.min_loop_length)
        end = time.time()
        logging.debug(f'Razor done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('loops acquired with razor'))
        print("=" * 50)

    # 8.Autoimmunity/allergenicity module
    if args.autoimmunity == "True":
        print("=" * 50)
        print("{:^50}".format('Autoimmunity and allergenicity starts'))
        print("=" * 50)
        logging.debug("Autoimmunity and allergen module starts ...")
        start = time.time()

        list_of_proteins = autoimmunity(list_of_proteins, args.proteome1, args.working_dir, args.NERVE_dir, args.e_value,
                                        args.minlength, args.mismatch, args.substitution)

        end = time.time()
        logging.debug(f'human Autoimmunity done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('Human autoimmunity and allergenicity acquired'))
        print("=" * 50)

    # mouse immunity and allergenicity module
        if args.mouse == "True":
            print("=" * 50)
            print("{:^50}".format('Mouse autoimmunity and allergenicity starts'))
            print("=" * 50)
            start = time.time()
            logging.debug('Mouse immunity starts ...')
            list_of_proteins = mouse(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value, args.proteome1,
                                 args.minlength, args.substitution, args.mismatch)
            end = time.time()
            logging.debug(f'mouse autoimmunity Done in {end - start} seconds')
            print("=" * 50)
            print("{:^50}".format('Mouse autoimmunity and allergenicity acquired'))
            print("=" * 50)

    # 9.Conservation pr1/pr2 module
    if args.proteome2:
        print("=" * 50)
        print("{:^50}".format('Conservation starts'))
        print("=" * 50)
        start = time.time()
        logging.debug("Conservation module starts ...")
        list_of_proteins = conservation(list_of_proteins, args.working_dir, args.NERVE_dir, args.e_value,
                                        args.proteome1, args.proteome2, args.minlength, args.substitution,
                                        args.mismatch)
        end = time.time()
        logging.debug(f'Conservation done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('Conservation acquired'))
        print("=" * 50)

    # 10.Annotation module (deepFRI)
    if args.annotation == "True":
        print("=" * 50)
        print("{:^50}".format('Annotation with deepFRI starts'))
        print("=" * 50)
        start = time.time()
        logging.debug('Annotation module starts ...')
        list_of_proteins = annotation(list_of_proteins, args.proteome1, args.working_dir, args.DeepFri_dir)
        end = time.time()
        logging.debug(f'Annotation done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('Annotation acquired'))
        print("=" * 50)
        
    # calculate instability index and charge 
    list_of_proteins = protparam(list_of_proteins)

    # 11.Protein selection module
    final_proteins = list_of_proteins
    if args.select == "True":
        print("=" * 50)
        print("{:^50}".format('Protein selection starts'))
        print("=" * 50)
        logging.debug('Selection module starts ...')
        start = time.time()
        final_proteins = select(list_of_proteins, args.autoimmunity, args.transmem_doms_limit, args.padlimit, args.mouse,
                                args.antigenlimit, args.antigen, args.annotation, args.razor)
        end = time.time()
        logging.debug(f'Selection done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('Protein selection and scoring done'))
        print("=" * 50)

    # 12.Epitope prediction
    if args.epitopes == "True":
        print("=" * 50)
        print("{:^50}".format('Epitope prediction of best candidates with epitopepredict starts'))
        print("=" * 50)
        start = time.time()
        logging.debug('Epitope prediction starts ...')
        final_proteins = epitope(final_proteins, args.autoimmunity, args.mouse,
                                 args.antigen, args.working_dir, args.mhci_length, args.mhcii_length,
                                 args.mhci_overlap, args.mhcii_overlap, args.epitope_percentile)
        end = time.time()
        logging.debug(f'Epitope prediction done in {end - start} seconds')
        print("=" * 50)
        print("{:^50}".format('Epitopes acquired'))
        print("=" * 50)

    print("=" * 50)
    print("{:^50}".format('CSV generation starts'))
    print("=" * 50)

    # Return csv file with outputs
    output(final_proteins, os.path.join(args.working_dir, 'eu_vaccine_candidates.csv'),
           args.mouse, args.autoimmunity, args.antigen)
    # collect also discarded proteins
    final_proteins_names = [p.id for p in final_proteins]
    discarded_proteins = [p for p in list_of_proteins if p.id not in final_proteins_names]
    output(discarded_proteins, os.path.join(args.working_dir, 'discarded_proteins.csv'),
           args.mouse, args.autoimmunity, args.antigen)

    enerve_end = time.time()
    logging.debug(f'eNERVE has finished its analysis in: {enerve_end - enerve_start} seconds')
    print("=" * 50)
    print("{:^50}".format(f'CSV of candidate and discarded proteins saved in {args.working_dir}'))
    print("=" * 50)

    print("=" * 50)
    print("{:^50}".format(f'End of eNERVE computation'))
    print("=" * 50)


if __name__ == '__main__':
  main()

