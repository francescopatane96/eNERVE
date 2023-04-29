###
### eNERVE dependencies installation
###
### base image found here https://hub.docker.com/_/python
### first stage (base, with Python 3.7)


ARG VERSION
FROM python:3.7 AS base
 # set where NERVE is contained
ENV PYTHONPATH "${PYTHONPATH}:/eNERVE"
# set where spaan is contained
ENV PYTHONPATH "${PYTHONPATH}:/"
ENV LANG C

FROM base AS intermediate

FROM ubuntu:14.04.4
RUN apt-get update && apt-get install -y apt-transport-https
RUN pip install setuptools==65.4.1
RUN \
  apt-get update && \
  apt-get install -y sudo curl git && \
  curl -s https://packagecloud.io/install/repositories/github/git-lfs/script.deb.sh | sudo bash && \
  sudo apt-get install git-lfs=1.0.0
RUN git lfs install 
RUN git lfs clone https://github.com/FranceCosta/nerve_data.git

RUN mv /nerve_data/DeepFri.tar.gz DeepFri.tar.gz && \
    tar xvzf ./DeepFri.tar.gz

RUN mv /nerve_data/newest_trained_models.tar.gz trained_models.tar.gz && \
    tar xvzf trained_models.tar.gz -C ./DeepFri

RUN mv /nerve_data/iFeature.tar.gz iFeature.tar.gz && \
    tar xvzf iFeature.tar.gz

WORKDIR ./DeepFri
RUN pip install .
WORKDIR /

# install python dependencies
COPY requirements.txt .
RUN pip install -r ./requirements.txt && \
    python -m pip install git+https://github.com/nicolagulmini/tmhmm.py
RUN pip install tensorflow
RUN pip install -U scikit-learn

FROM intermediate AS dependencies

# clone repositories
#RUN git clone https://github.com/Superzchen/iFeature.git 
RUN  git clone https://github.com/francescopatane96/eNERVE
     

RUN apt-get update 
RUN apt-get install -y apt-utils ncbi-blast+ 
RUN apt-get install nano

FROM dependencies AS setting
ARG USER_ID
ARG GROUP_ID
EXPOSE 8880
RUN addgroup --gid $GROUP_ID user # name user as user
RUN adduser --disabled-password --gecos '' --uid $USER_ID --gid $GROUP_ID user # name group user as user
RUN mkdir /workdir
RUN chown -R user:user /workdir
# give permissions to allow testing
RUN chown -R user:user /eNERVE
RUN chmod +x /eNERVE/enerve/nerve.py
USER user

COPY . .
CMD ["bash"]

FROM setting as version-interactive
ENV cmd="jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''"

FROM setting as version-commandline
ENV cmd="['/eNERVE/code/nerve.py']"

#FROM version-${VERSION} AS final
#RUN echo "command is ${cmd}"
#CMD cmd
CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''

#CMD jupyter notebook --notebook-dir=/workdir --ip='0.0.0.0' --port='8880' --allow-root --NotebookApp.token='' --NotebookApp.password=''
#ENTRYPOINT ["/eNERVE/enerve/nerve.py"]

