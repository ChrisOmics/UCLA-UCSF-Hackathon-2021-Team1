FROM rocker/rstudio:latest

WORKDIR /home

COPY install_monocle_and_R_dependencies.R .

RUN Rscript install_monocle_and_R_dependencies.R

WORKDIR /data

COPY example-data .

WORKDIR /tutorial

COPY single-cell-tutorial.Rmd .

WORKDIR /home
