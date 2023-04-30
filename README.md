# eNERVE v1.0 - Eucaryotic New Enhanced Reverse Vaccinology Environment

 eNERVE is a dynamic, high-throughput and standalone in silico pipeline for PVCs (protein vaccine candidates) discovery in eukaryotic organisms, by reverse vaccinology methods. 
 The tool is capable of identifying candidate antigens (immunogens) and their epitopes with ML (machine learning) and alignment analysis methods, using a tree based approach, from entire proteomes as input (FASTA file of proteomes downloaded from Uniprot and other databases or generated from proteome experiments).
 The output generated by the pipeline is a general CSV file containing every candidate protein with relative score and predicted values (P_adhesin, P_antigen, location etc..), a CSV file with discarded proteins and a series of files (CSVs, PNGs) located in specific folders for every protein selected to predict coresponding epitopes. 
 
 eNERVE is a flexible tool and permits you to select specific modules and cutoffs to use during proteome analysis and protein classification task (please, see ```Usage``` section of this README for more informations).

eNERVE is designed to assist experimental research activities in vaccine discovery and vaccine formulation for eucaryotic targets, making use if the data available to the scientific community and using it to create machine learning models that can facilitate and economize the process of vaccine antigen discovery.

Are you searching for a bacteria vaccine discovery pipeline? Please visit https://nicolagulmini.github.io/NERVE/ and https://github.com/FranceCosta/NERVE .
 
 ***
 # Pipeline architecture and data flow:
 
 1. Quality control of proteome module and generation of protein instances: proteins that do not pass QC process are discarded (discarded_sequences.fasta);
 2. Descriptors calculator module: this module uses iFeature library (https://github.com/Superzchen/iFeature) to calculate protein descriptors for every protein in the input proteome;
 3. Subcellular location module: RandomForest based model (scikit-learn) to predict the probability of being 'outer' (training dataset from UniProt https://www.uniprot.org/);
 4. Transmembrane topology prediction module: TMHMM (Hidden markov model based model), predicts the probability of each aminoacid to be 'i', 'o' or 'M'. For more informations, please visit https://github.com/dansondergaard/tmhmm.py ;
 5. Razor module: virtual scissors that the outer protein pieces and rejoin them. This is useful to retrieve the outer segments of proteins with many transmembrane domains and to discard the latter.
 6. Adhesin and adhesin-like predictor module: feed-forward neural network (Tensorflow/Keras). Training dataset obtained from literature and InterPro (https://www.ebi.ac.uk/interpro/);
 7. Antigenicity predictor module: feed-forward neural network (Tensorflow/Keras). Training data obtained from IEDB database (https://www.iedb.org/)
 8. Autoimmunity and allergenicity module: Alignment based method (ncbi blast+)
 9. Selection and scoring module
 10. Linear epitope predictor module: epitopepredict (https://github.com/dmnfarrell/epitopepredict)
 11. Conservation module between proteome1 and proteome2, if proteome2 is added to the analysis
 12. Output generation module
 ***
 
 ### Instructions for stand-alone usage with docker: under development
eNERVE can be used as a stand-alone version taking advantage of [Docker](https://www.docker.com/) and [Docker-compose](https://docs.docker.com/engine/reference/commandline/compose/) in linux systems.

1) install Docker following [these instructions](https://docs.docker.com/engine/install/) and [the post-installation procedure](https://docs.docker.com/engine/install/linux-postinstall/)
2) install docker-compose as explained [here](https://docs.docker.com/compose/install/linux/)
3) clone repository:
```
git clone git@github.com:/francescopatane96/eNERVE.git
```
4) cd to docker folder
```
cd ./eNERVE/docker
```
5) build docker containers. This takes a few minutes
```
./build_run.sh
```
6) open a browser and navigate to local host: http://localhost:8880

 
 ***
 # Instructions for local installation (tested with python3.10.6 and 3.7.14 versions in linux/ubuntu/unix os environment):
 
 1. Open the terminal, move to the location in which you would save eNERVE (eg. ```cd /home/ubuntu/Desktop/```);
 2. ``` git clone https://github.com/francescopatane96/eNERVE.git ``` ;
 3. ``` cd eNERVE/code ```;
 4. ``` git clone https://github.com/francescopatane96/DeepFRI.git ```;
 5. Download pre-trained models from https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz, then uncompress ```tar.gz file``` into the ```DeepFRI directory``` (```tar xvzf trained_models.tar.gz -C /path/to/DeepFRI```);
 6. ``` git clone https://github.com/francescopatane96/iFeature.git ``` and ```pip install numpy```;
 7. ``` pip install git+https://github.com/francescopatane96/tmhmm.py.git ```;
 8. ``` sudo apt-get install ncbi-blast+ ```;
 9. ``` pip install tensorflow ```;
 10. ``` pip install -U scikit-learn ```;
 11. finally, ``` pip install -r requirements.txt ```
 
***
# Instructions for creating a virtual environment (linux/unix/mac os):
1. open a terminal or a terminal from an IDE (pyCharm or visual code studio)
2. ```cd /path/destination/``` and clone the repository (```git clone ...```)
3. ```cd eNERVE```
4. Install ```python3-venv``` package digiting ```sudo apt install python3.10-venv ```
5. create a virtual environment with python module ```venv``` (to avoid dependencies conflicts) with ```python3 -m venv enerve```
6. activate your new virtual environment with ```source enerve/bin/activate```
After venv activation, the terminal will shows virtual environment name between ```()```, eg. ```(nerve)```
7. ```cd eNERVE/code```
8. Now, you have to install dependencies (from point 4 of the previous section) needed for the pipeline.

***
# Usage:
```
usage: nerve.py [-h] [-a] [-ai] [-tp] [-ev] [-ml] [-mm] [-m]
                [-mpsl] -p1 [-p2] [-pl] [-rz] [-ig] [-rl] [-s]
                [-ss] [-tdl] [-ang] [-wd] [-nd] [-id] [-dfd]
                [-ep] [-m1l] [-m2l] [-m1ovr] [-m2ovr] [-prt]
                
                where:
                -h (help), [];
                -a (Protein functional annotation with DeepFRI), [True, False, default=True];
                -ai (autoimmunity and allergenicity module), [True, False, default=True];
                -tp (topology), [tmhmm];
                -ev (e-value for blastp), [float, default=1e-10];
                -ml (minlength required for shared peptides to be extracted in comparison analysis versus human and/or mouse) [int, default=9];
                -mm (mismatch, maximal number of not compatible substitutions allowed in shared peptides alignment windows of minlength size in immunity module, [int, default=1];
                -m (mouse autoimmunity and allergenicity module), [True, False, default=True];
                -mpsl (mouse peptides sum limit, parameter used by selection module. protein with sum of shared peptides of the i protein with mouse proteins/number of aminoacids of the i protein <= mouse peptides sum limit and with absence of match mhc-I and mhc-II mouse ligands are selected), [float, default=0.15];
                -p1 (proteome 1 fasta filename or path), [filename.fasta] --> required;
                -p2 (proteome 2 fasta filename or path), [filename.fasta];
                -rz (razor module), [True, False, default=True];
                -ig (antigenlimit, cutoff value for antigen module), [float, default=0.80];
                -rl (min loop length considered in razor module), [int, default=9];
                -s (selection module), [True, False, default=True];
                -ss (substitution, maximal number of compatible substitutions allowed in shared peptides alignment windows of minlength size in immunity module), [int, default=3];
                -tdl (transmembrane doms limit) [int, default=30];
                -ang (antigen module), [True, False, default=True];
                -wd [path/to/workdir] --> recommended;
                -nd [path/to/NERVEdir] --> recommended;
                -id (iFeature directory), [path/to/ifeature_dir, default=./iFeature];
                -dfd (DeepFri directory), [path/to/deepfri_dir, default=./DeepFRI];
                -ep (epitope prediction module), [True, False, default=True];
                -m1l (mhc1 ligands length), [9,10,11, default=9];
                -m2l (mhc2 ligands length), [9,11,13,15, default=11];
                -m1ovr (mhc1 ligands max overlap), [1,2,default=1];
                -m2ovr (mhc2 ligands max overlap), [1,2, default=1];
                -prt (epitope binders percentile), [float, default=0.9]
```

 1. digit ```cd eNERVE```;
 2. create your working directory (where you will put in fasta.file to be analyze and where outputs will be saved);
 3. ```cd code```;
 3. In the terminal, digit and run ```python3 nerve.py``` followed by args**;
 4. Output files will be saved in the working directory at the end of the computation
 
 
 ***
# References and contacts:
 eNERVE was developed by Francesco Patanè during his master thesis and internship under the supervision of Prof. Francesco Filippini, at University of Padova.
 
 Special thanks to Francesco Costa, Nicola Gulmini and Andrea Conte for their help and collaboration in this project.
 
 This pipeline was also implemented through the use of packages and libraries created by others (iFeature, epitopepredict, tmhmm, ncbi-blast+, tensorflow and many others), so thanks to https://github.com/dansondergaard, https://github.com/dmnfarrell/epitopepredict, https://github.com/francescopatane96/iFeature and thanks to all the open source community for their pretty and useful tools to predict transmembrane protein topology, linear epitopes and so on.
 
 
