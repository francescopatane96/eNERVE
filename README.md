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
 3. Subcellular location module: RandomForest based model to predict the probability of being 'outer' (training dataset from UniProt https://www.uniprot.org/);
 4. Transmembrane topology prediction module: TMHMM (Hidden markov model based model), predicts the probability of each aminoacid to be 'i', 'o' or 'M'. For more informations, please visit https://github.com/dansondergaard/tmhmm.py ;
 5. Razor module: virtual scissors that the outer protein pieces and rejoin them. This is useful to retrieve the outer segments of proteins with many transmembrane domains and to discard the latter.
 6. Adhesin and adhesin-like predictor module: feed-forward neural network (Tensorflow/Keras). Training dataset obtained from literature and InterPro (https://www.ebi.ac.uk/interpro/);
 7. Antigenicity predictor module: feed-forward neural network (Tensorflow/Keras)
 8. Autoimmunity and allergenicity module: Alignment based method
 9. Selection and scoring module
 10. Linear epitope predictor module: Epitopepredict
 11. Conservation module between proteome1 and proteome2, if proteome2 is added to the analysis
 12. Output generation module
 ***
 
 # Instructions for stand-alone installation with docker: under development
 
 ***
 # Instructions for local installation (tested with python3.10 and 3.7 versions in linux/ubuntu/unix os environment):
 
 1. Open the terminal, move to the location in which you would save eNERVE (eg. ```cd /home/ubuntu/Desktop/```);
 2. ``` git clone https://github.com/francescopatane96/eNERVE.git ``` ;
 2. ``` cd eNERVE/code ```;
 3. ``` git clone https://github.com/francescopatane96/DeepFRI.git ```;
 3. Download pre-trained models from https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz, then uncompress ```tar.gz file``` into the ```DeepFRI directory``` (```tar xvzf trained_models.tar.gz -C /path/to/DeepFRI```);
 4. ``` git clone https://github.com/francescopatane96/iFeature.git ``` and ```pip install numpy```;
 5. ``` pip install git+https://github.com/francescopatane96/tmhmm.py.git ```;
 6. ``` sudo apt-get install ncbi-blast+ ```;
 7. ``` pip install tensorflow ```;
 8. ``` pip install -U scikit-learn ```;
 9. finally, ``` pip install -r requirements.txt ```
***
# Instructions for creating a virtual environment (linux/unix/mac os):
1. open a terminal or a terminal from an IDE (pyCharm or visual code studio)
2. ```cd /path/destination/``` and clone the repository (```git clone ...```)
3. ```cd eNERVE/code```
4. create a virtual environment with python module ```venv``` (to avoid dependencies conflicts) with ```python3 -m venv enerve```
5. activate your new virtual environment with ```source enerve/bin/activate```
After venv activation, the terminal will shows virtual environment name between ```()```, eg. ```(nerve)```
6. Now, you have to install dependencies (from point 2 of the previous section) needed for the pipeline.

***
# Usage:
 1. In the terminal, digit and run ```python3 nerve.py -arg1 -arg2 -arg(n+1)```
 
 
 ***
# References and contacts:
 eNERVE was developed by Francesco Patanè during his master thesis and internship under the supervision of Prof. Francesco Filippini, at University of Padova.
 Thanks to Francesco Costa, Nicola Gulmini and Andrea Conte for their help and collaboration in this project, and to https://github.com/dansondergaard and https://github.com/dmnfarrell/epitopepredict for their pretty and useful tools to predict transmembrane protein topology and linear epitopes.
 
 
