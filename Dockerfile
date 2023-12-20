FROM condaforge/mambaforge:latest

ADD environment.yaml /

RUN /opt/conda/bin/mamba env create --name masq --file /environment.yaml

# GPF ENV
ENV PATH /opt/conda/envs/masq/bin:$PATH

RUN mkdir -p /data && mkdir -p /code && mkdir -p /wd

WORKDIR /wd

SHELL ["/bin/bash", "-c"]
