# eNERVE
 A dynamic, high-throughput and standalone in silico pipeline for eukaryotic reverse vaccinology. 
 eNERVE is a tool for reverse vaccinology capables of performing ML (machine learning) and alignment analysis over a tree based approach from entire proteomes as input (FASTA file).
 
 1. Quality control module
 2. Sublocation module: RandomForest based model
 3. Topology module: TMHMM
 4. Adhesin and adhesin-like predictor module: feed-forward neural network
 5. Antigenicity module
 6. Autoimmunity and allergenicity module: Alignment based
 7. Selection and scoring module
 8. Linear epitope predictor module
 9. Conservation module
 
 
 Instructions for expert users:
 
 1. open the terminal, change and open the location in which you would save eNERVE.
 2. git clone https://github.com/francescopatane96/eNERVE.git
 2. cd eNERVE
 3. git clone https://github.com/francescopatane96/DeepFRI.git
 4. git clone https://github.com/francescopatane96/iFeature.git
 5. pip install git+https://github.com/francescopatane96/tmhmm.py.git
 6. apt-get install ncbi-blast+
 7. pip install tensorflow
 8. pip/pip3 install -U scikit-learn
 9. finally, pip install -r requirements.txt
 
 Usage:
 
