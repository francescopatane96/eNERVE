# eNERVE v1.0 - Eucaryotic New Enhanced Reverse Vaccinology Environment

[![run with conda](http://img.shields.io/badge/run%20with-conda-3EB049?labelColor=000000&logo=anaconda)](https://docs.conda.io/en/latest/)
[![run with docker](https://img.shields.io/badge/run%20with-docker-0db7ed?labelColor=000000&logo=docker)](https://www.docker.com/)
![Jupyter Notebook](https://img.shields.io/badge/jupyter-%23FA0F00.svg?style=for-the-badge&logo=jupyter&logoColor=white)
![Python](https://img.shields.io/badge/python-3670A0?style=for-the-badge&logo=python&logoColor=ffdd54)
![Keras](https://img.shields.io/badge/Keras-%23D00000.svg?style=for-the-badge&logo=Keras&logoColor=white)
![Pandas](https://img.shields.io/badge/pandas-%23150458.svg?style=for-the-badge&logo=pandas&logoColor=white)
![NumPy](https://img.shields.io/badge/numpy-%23013243.svg?style=for-the-badge&logo=numpy&logoColor=white)
![Matplotlib](https://img.shields.io/badge/Matplotlib-%23ffffff.svg?style=for-the-badge&logo=Matplotlib&logoColor=black)
![TensorFlow](https://img.shields.io/badge/TensorFlow-%23FF6F00.svg?style=for-the-badge&logo=TensorFlow&logoColor=white)
![Linux](https://img.shields.io/badge/Linux-FCC624?style=for-the-badge&logo=linux&logoColor=black)
![macOS](https://img.shields.io/badge/mac%20os-000000?style=for-the-badge&logo=macos&logoColor=F0F0F0)
![Ubuntu](https://img.shields.io/badge/Ubuntu-E95420?style=for-the-badge&logo=ubuntu&logoColor=white)
![Windows](https://img.shields.io/badge/Windows-0078D6?style=for-the-badge&logo=windows&logoColor=white)
<img alt="Issues" src="https://camo.githubusercontent.com/3093fdb6f8038165918c5338bb94d4b9983114e635e34f396c03b3063ff206fc/68747470733a2f2f696d672e736869656c64732e696f2f6769746875622f6973737565732f77696e7465726c6f6f642f726561646d652d737469636b6572733f636f6c6f723d303038386666" data-canonical-src="https://img.shields.io/github/issues/winterlood/readme-stickers?color=0088ff" style="max-width: 100%;">



eNERVE is a dynamic, high-throughput and standalone in silico pipeline for ```PVCs``` (protein vaccine candidates) discovery in ```eukaryotic``` organisms, by a ```reverse vaccinology``` approach. 
 
The tool is capable of identifying candidate ```antigens``` (immunogens) and their ```epitopes``` with ```ML``` (machine learning) and ```alignment``` analysis methods, using a tree based approach, from entire proteomes as input (```FASTA``` file of proteomes downloaded from Uniprot and other databases or generated from proteome experiments).

The output generated by the pipeline is a general ```CSV``` file containing every candidate protein with relative score and predicted values (P_adhesin, P_antigen, location etc..), a CSV file with discarded proteins and a series of files (CSVs, PNGs) located in specific folders for every protein selected to predict corresponding epitopes. 
 
eNERVE is a flexible tool that allows you to select specific modules and cutoffs to use during proteome analysis and protein classification task (please, see ```Usage``` section of this ```README``` for more informations).

eNERVE is designed to assist experimental research activities in vaccine discovery and vaccine formulation for eucaryotic targets, making use of the data available to the scientific community and using it to create machine learning models that can facilitate and economize the process of antigens discovery to formulate protein and subunit vaccines.

Are you searching for a ```bacteria``` vaccine discovery pipeline? Please visit NERVE (```bNERVE```) repository, [here](https://nicolagulmini.github.io/NERVE/) and [here]( https://github.com/FranceCosta/NERVE) .
 
 ***
 ### 💻 Pipeline architecture and data flow:
 ![alt text](https://github.com/francescopatane96/eNERVE/blob/main/workflow_design.png)

 1. Proteins or entire proteomes downloaded by Databases and proteomic experiments and in .FASTA format are passed to the pipeline;
 2. Quality control of proteome module and generation of protein instances: proteins that do not pass QC process are discarded (discarded_sequences.fasta);
 Descriptors are calculated using [iFeature library](https://github.com/Superzchen/iFeature) to generate protein descriptors for every protein in the input proteome;
 3. Subcellular location module: Random Forest-based model (scikit-learn) to predict the probability of being 'outer' (training dataset from [UniProt]( https://www.uniprot.org/));
 4. Adhesin and adhesin-like predictor module: feed-forward neural network (Tensorflow/Keras). Training dataset obtained from literature and [InterPro](https://www.ebi.ac.uk/interpro/);
 5. Autoimmunity and allergenicity module: Alignment based method (ncbi blast+). Input proteins are aligned with human and mouse proteome in order to catch the level of similarity. What is more, a list of autoimmune and allergenic peptides are screened on the input proteome;
 6. Transmembrane topology prediction module: TMHMM (Hidden markov model). Predicts the probability of each aminoacid to be 'i', 'o' or 'M'. For more informations, please visit tmhmm [repository](https://github.com/dansondergaard/tmhmm.py) ;
 7. Razor module: virtual scissors that cut outer protein pieces and rejoin them. This is useful to retrieve the outer segments of proteins with many transmembrane domains and to discard the latter.
 9. Conservation module: if also the proteome2 is added to the analysis, proteomes are compared.
 10. Selection and scoring module: external proteins are saved in the "vaccine_candidates" file, while internal proteins are saved in the "discarded_proteins" one. Internal proteins with a P_ad > padlimit threshold are retained and saved in the candidates file.
 11. Linear epitope predictor module: [epitopepredict library](https://github.com/dmnfarrell/epitopepredict). For every protein the module predicts its linear epitopes and promiscuous epitopes considering only the 'supertypes alleles' defined by Sette et al.;
 12. Output generation module: a .CSV file is generated which contains every protein instance with all predictions (columns) like score, p_ad, p_loc_out, transmem doms, epitopes, lenght, instability index etc...
 ***
 
 ### :accessibility: Instructions for stand-alone usage with Docker and dockerhub (preferred method):

eNERVE can be used as a stand-alone version taking advantage of [Docker](https://www.docker.com/) and [Docker-compose](https://docs.docker.com/engine/reference/commandline/compose/) in linux systems. This method ensures no dependencies related issues.

1) install Docker following [these instructions](https://docs.docker.com/engine/install/) and [the post-installation procedure](https://docs.docker.com/engine/install/linux-postinstall/)
2) install docker-compose as explained [here](https://docs.docker.com/compose/install/linux/)
3) open the terminal and digit:
```
sudo docker pull francescopatane/enerve:tag
```
4) create a directory (eg. on Desktop) called 'output' and give it permissions with:
```
chmod 777 /path/to/output
```
5) Then, put your input FASTA files in the directory;

6) Run docker image and select a volume for sharing input from local machine and output from virtual machine:
```
sudo docker run --rm -it -v /path/to/output_directory:/workdir francescopatane/enerve:tag
```
7) go to root folder with:
```
cd ..
```
8) run eNERVE pipeline with:
```
python3 nerve.py -wd [], -dfd [./], -p1 [filename.fasta] -args**
```
9) At the end of the computation you will find output files in your 'output' directory in the local machine, in this case in your Desktop in the directory ```output```

 
 ***
 ### 🏠 Instructions for local installation (tested with python3.10.6 and 3.7.14 versions in Linux/Ubuntu/Unix OS environments):
 
 1. Open the terminal, move to the location in which you would save eNERVE:
 ```
 cd /home/ubuntu/Desktop/
 ```
 2. clone the repository on your machine:
 ``` 
 git clone https://github.com/francescopatane96/eNERVE.git
 ``` 
 3. move to the directory:
 ``` 
 cd eNERVE
 ```
 4. Clone DeepFri repository:
 ``` 
 git clone https://github.com/francescopatane96/DeepFRI.git 
 ```
 5. Download pre-trained models from [flatironinstitute](https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz), then uncompress ```tar.gz file``` into the ```DeepFRI directory```:
 ```
 tar xvzf newest_trained_models.tar.gz -C /path/to/DeepFRI
 ```
 6. clone iFeature repository:
 ``` 
 git clone https://github.com/francescopatane96/iFeature.git
 ``` 
 and install numpy:
 ```
 pip install numpy
 ```
 7. install tmhmm library:
 ``` 
 pip install git+https://github.com/francescopatane96/tmhmm.py.git 
 ```
 8. install ncbi-blast +:
 ``` 
 sudo apt-get install ncbi-blast+ 
 ```
 9. install tensorflow:
 ``` 
 pip install tensorflow 
 ```
 10. install scikit-learn library:
 ``` 
 pip install -U scikit-learn 
 ```
 11. finally, install remaining dependencies:
 ``` 
 pip install -r requirements.txt 
 ```
 
***
### 🤖 Instructions for creating a virtual environment with Python venv (Linux/Unix/Mac OS):
1. open a terminal or a terminal from an IDE (pyCharm or visual code studio);
2. move to destination folder:
```
cd /path/destination/
``` 
and clone the repository :
```
git clone https://github.com/francescopatane96/eNERVE.git
```
3. move to root directory:
```
cd eNERVE
```
4. Install ```python3-venv``` package digiting:
```
sudo apt install python3.10-venv
```
5. create a virtual environment with python module ```venv``` (to avoid dependencies conflicts) with: 
```
python3 -m venv enerve
```
6. activate your new virtual environment with: 
```
source enerve/bin/activate
```
After venv activation, the terminal will shows virtual environment name between ```()```, eg. ```(nerve)```;

7. Now, you have to install dependencies (from point 4 of the previous section) needed for the pipeline.
***
### 🔴 We recommend using conda for creating a virtual environment:

1. from the terminal, digit: 
```
conda create --name enerve python=3.10
```
2. Clone the repository;

3. activate the environment:
```
conda activate enerve
```
4. install dependencies as in previous sections.

***
### 🦮 Usage:
```
usage: nerve.py [-h] [-locl] [pl] [-a] [-ai] [-tp] [-ev] [-ml] [-mm] [-m]
                [-mpsl] -p1 [-p2] [-rz] [-ig] [-rl] [-s]
                [-ss] [-tdl] [-ang] [-wd] [-nd] [-id] [-dfd]
                [-ep] [-m1l] [-m2l] [-m1ovr] [-m2ovr] [-prt]
                
                where:
                -h (help), [];
                -locl (loclimit, localization prediction threshold outer class), [float, default=0.60],
                -pl (adhlimit, Retrieve internal proteins if having adh probability > adl), [float, default=0.80];
                -a (Protein functional annotation with DeepFRI), [True, False, default=True];
                -ai (human autoimmunity and allergenicity module), [True, False, default=True];
                -tp (topology), [tmhmm];
                -ev (e-value for blastp), [float, default=1e-10];
                -ml (minlength required for shared peptides to be extracted in comparison analysis versus human and/or mouse) [int, default=9];
                -mm (mismatch, maximal number of not compatible substitutions allowed in shared peptides alignment windows of minlength size in immunity module, [int, default=1];
                -m (mouse autoimmunity and allergenicity module), [True, False, default=True];
                -mpsl (mouse peptides sum limit, parameter used by selection module. protein with sum of shared peptides of the i protein with mouse proteins/number of aminoacids of the i protein <= mouse peptides sum limit and with absence of match mhc-I and mhc-II mouse ligands are selected), [float, default=0.15];
                -p1 (proteome 1 fasta filename or path), [filename.fasta] --> 🔴required🔴;
                -p2 (proteome 2 fasta filename or path), [filename.fasta];
                -rz (razor module), [True, False, default=True];
                -ig (antigenlimit, cutoff value for antigen module), [float, default=0.80] --> 🔵Experimental/Beta🔵;
                -rl (min loop length considered in razor module), [int, default=9];
                -s (selection module), [True, False, default=True];
                -ss (substitution, maximal number of compatible substitutions allowed in shared peptides alignment windows of minlength size in immunity module), [int, default=3];
                -tdl (transmembrane doms limit) [int, default=0] --> 🔵For whole protein vaccines, use a number != 0. For epitope-based predictions, use 0🔵;
                -ang (antigen module), [True, False, default=False] --> 🔵Experimental/Beta🔵;
                -wd [path/to/workdir] --> 🟠recommended🟠;
                -nd [path/to/NERVEdir] --> 🟠recommended🟠;
                -id (iFeature directory), [path/to/ifeature_dir, default=./iFeature];
                -dfd (DeepFri directory), [path/to/deepfri_dir, default=./DeepFRI];
                -ep (epitope prediction module), [True, False, default=True];
                -m1l (mhc1 ligands length), [9,10,11, default=9];
                -m2l (mhc2 ligands length), [9,11,13,15, default=11];
                -m1ovr (mhc1 ligands max overlap), [1,2,default=1];
                -m2ovr (mhc2 ligands max overlap), [1,2, default=1];
                -prt (epitope binders percentile), [float, default=0.80]

```
:warning:
Remember that required and essential parameters are [-wd], [-p1] and [-nd] if using Local installation. Dockerized version needs only [-wd]. By default, every module is active and will be run. To personalized and deactivate single modules, digit -parameter** False 
:warning:
***
### ⛑️ Examples

# To run eNERVE in a local environment or in a virtual one (no docker) with all modules:
Digit on the command line:

```
python3 code/Nerve.py -p1 proteome.fasta
```
To run eNERVE without annotation (-a) module, mhc1 ligands length of 10, mhc2 ligands length of 15 and epitope percentile of 80:
```
python3 code/Nerve.py -p1 proteome.fasta -a False -prt 0.80 -m1l 10 -m2l 15
```
If you want to create a workid directory in which save outputs, please specify [-wd] (eg. -wd /workdir)

 1. digit:
 ```
 cd eNERVE
 ```
 2. create your working directory (where you will put in fasta.file to be analyze and where outputs will be saved);
 
 3. In the terminal, digit and run:
 ```
 python3 code/Nerve.py -arg1 -arg2 -args**
 ```
:warning:
 REMEMBER TO SPECIFY -wd (WORKING_DIR) if desired, -p1 (eventually also -p2). Place your fasta inputs into the workdir if just exist.
 :warning:
 
 4. Output files will be saved in the working directory at the end of the computation

# If you are using Docker version of eNERVE:

 
 
 ***
### 📲 References and contacts:
 eNERVE was developed by ```Francesco Patanè``` during his master thesis and internship under the supervision of Prof. ```Francesco Filippini```, at ```University of Padova```, Synthetic Biology and Biotechnology unit (SynBio) .
 
 Special thanks to [Francesco Costa](https://github.com/FranceCosta), [Nicola Gulmini](https://github.com/nicolagulmini) and ```Andrea Conte``` for their help and collaboration in this project.
 
 This pipeline was also implemented through the use of packages and libraries created by others (iFeature, epitopepredict, tmhmm, ncbi-blast+, tensorflow and many others), so thanks to [dansondergaard](https://github.com/dansondergaard), [dmnfarrell](https://github.com/dmnfarrell/epitopepredict), [superzchen](https://github.com/Superzchen/iFeature) for their pretty and useful tools to predict transmembrane protein topology, linear epitopes and thanks to all the open source community (in particular, the machine learning one).
 
Have you encountered any problems installing or using the pipeline, or have any suggestions for improving eNERVE? please contact me at the following addresses:

 ```
 francesco.patane@live.it
 francesco.patane.1@studenti.unipd.it
 ```
 or open an [issue](https://github.com/francescopatane96/eNERVE/issues) 

 
 ```
1. Vivona S, Bernante F, Filippini F. NERVE: new enhanced reverse vaccinology environment. BMC Biotechnol. 2006 Jul 18;6:35. doi: 10.1186/1472-6750-6-35. PMID: 16848907; PMCID: PMC1570458.
 ```
<!--

Thesis project by [MOLBINFO](http://www.bio.unipd.it/molbinfo/). 
[eNERVE](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC1570458/) means '*New Enhanced Reverse Vaccinology Environment*', the aim of the eNERVE project is to develop a high performance pipeline for selecting best candidate vaccines to target eukaryotic pathogenic organisms and to update single modules and training sets with new data available in the literature and from public databases.

-->

 
