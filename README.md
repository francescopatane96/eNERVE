# eNERVE v1.0 - Eucaryotic New Enhanced Reverse Vaccinology Environment

 eNERVE is a dynamic, high-throughput and standalone in silico pipeline for PVCs (protein vaccine candidates) discovery in eukaryotic organisms, by reverse vaccinology methods. 
 The tool is capable of identifying candidate antigens with ML (machine learning) and alignment analysis methods, using a tree based approach, from entire proteomes as input (FASTA file).
 
 1. Quality control of proteome module
 2. Subcellular location module: RandomForest based model
 3. Transmembrane topology prediction module: TMHMM
 4. Razor module
 5. Adhesin and adhesin-like predictor module: feed-forward neural network (Tensorflow/Keras)
 6. Antigenicity predictor module: feed-forward neural network (Tensorflow/Keras)
 7. Autoimmunity and allergenicity module: Alignment based method
 8. Selection and scoring module
 9. Linear epitope predictor module: Epitopepredict
 10. Conservation module
 
 
 Instructions for expert users (tested with python3.10 and 3.7 versions in linux/ubuntu/unix os environment)
 
 1. Open the terminal, change and open the location in which you would save eNERVE (eg. ```cd /home/ubuntu/Desktop/```);
 2. ``` git clone https://github.com/francescopatane96/eNERVE.git ``` ;
 2. ``` cd eNERVE/code ```;
 3. ``` git clone https://github.com/francescopatane96/DeepFRI.git ```;
 3. Download pre-trained models from ```https://users.flatironinstitute.org/~renfrew/DeepFRI_data/newest_trained_models.tar.gz```, then uncompress ```tar.gz file``` into the ```DeepFRI directory``` (```tar xvzf trained_models.tar.gz -C /path/to/DeepFRI```);
 4. ``` git clone https://github.com/francescopatane96/iFeature.git ```;
 5. ``` pip install git+https://github.com/francescopatane96/tmhmm.py.git ```;
 6. ``` apt-get install ncbi-blast+ ```;
 7. ``` pip install tensorflow ```;
 8. ``` pip install -U scikit-learn ```;
 9. finally, ``` pip install -r requirements.txt ```
 
 Usage:
 1. In the terminal, digit ```python3 nerve.py -arg1 -arg2 -arg(n+1) ... ```
 
 
 
 References and contacts:
 eNERVE was developed by Francesco Patanè during his master thesis and internship under the supervision of Prof. Francesco Filippini, at University of Padova.
 Thanks to Francesco Costa, Nicola Gulmini and Andrea Conte for their help and collaboration in this project.
 
 
