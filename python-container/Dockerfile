FROM continuumio/miniconda3:latest

WORKDIR /home

COPY env.yaml .


RUN conda update -n base -c defaults conda &&\
    conda env create -f env.yaml -n single-cell &&\
    echo "conda activate single-cell" >> ~/.bashrc &&\
    conda clean -afy &&\
    rm env.yaml

ENV PATH=/opt/conda/envs/single-cell/bin/:${PATH}

WORKDIR /home/tutorial

COPY example-data .
COPY python-tutorial.ipynb .

WORKDIR /home

CMD ["jupyter", "lab", "--allow-root", "--port=8888", "--no-browser", "--ip=0.0.0.0"]
