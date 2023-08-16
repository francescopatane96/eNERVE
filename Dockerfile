###
### eNERVE dependencies installation
###
### base image found here https://hub.docker.com/_/python
### first stage (base, with Python 3.10)


ARG VERSION
FROM python:3.10 AS base
 # set where NERVE is contained
ENV PYTHONPATH "${PYTHONPATH}:/eNERVE"

ENV LANG C

FROM base AS intermediate

#FROM ubuntu:14.04.4
RUN apt-get update && apt-get install -y apt-transport-https
RUN pip install setuptools==65.4.1
RUN apt-get install -y apt-transport-https
RUN apt-get update
RUN apt-get install -y sudo curl git 
RUN curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash
RUN sudo apt-get install git-lfs
RUN sudo git lfs install 
RUN sudo git lfs clone https://github.com/francescopatane96/nerve_data.git

WORKDIR ./DeepFri

RUN mv /nerve_data/DeepFri.tar.gz DeepFri.tar.gz && \
    tar xvzf ./DeepFri.tar.gz

RUN mv /nerve_data/newest_trained_models.tar.gz trained_models.tar.gz && \
    tar xvzf trained_models.tar.gz
RUN rm -r trained_models.tar.gz
RUN mv trained_models/* DeepFri && rm -r trained_models

RUN mv /nerve_data/iFeature.tar.gz iFeature.tar.gz && \
    tar xvzf iFeature.tar.gz
    
RUN sudo apt-get update
RUN pip install --upgrade pip


WORKDIR /

# install python dependencies
COPY requirements.txt .
RUN pip install requests 
RUN pip install -r ./requirements.txt && \
    python -m pip install git+https://github.com/nicolagulmini/tmhmm.py


FROM intermediate AS dependencies

# clone repositories
 
RUN git clone https://github.com/francescopatane96/eNERVE
RUN git clone https://github.com/francescopatane96/enerve_code.git
RUN mv enerve_code/* . && rm -r enerve_code
RUN mv eNERVE/* . && rm -r eNERVE
	
RUN apt-get update 
RUN apt-get install -y apt-utils ncbi-blast+ 
RUN apt-get install nano

FROM dependencies AS setting


# Create workingdir
WORKDIR /workdir

# Create a volume to share files
VOLUME ["/workdir"]

# Create new user
RUN useradd -ms /bin/bash newuser
RUN chown -R newuser:newuser /workdir
# change user
USER newuser

RUN chmod 755 /workdir


CMD ["bash"]




