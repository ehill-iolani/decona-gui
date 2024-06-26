# Set the base image; pretty sure rocker/shiny was built on Ubuntu 22.04 jammy
# not gamer, but fixing chip erra erra (like a dj), on arm devices like mac with m1/2 cheetos
FROM --platform=linux/amd64 rocker/shiny

# My authorship
LABEL maintainer="ehill@iolani.org"
LABEL version="1.0.0"
LABEL description="Decona GUI for Iolani School"

# Disable prompts during package installation
ENV DEBIAN_FRONTEND noninteractive

# Convenience packages
RUN apt update
RUN apt install -y curl git g++ zlib1g-dev make bsdmainutils gawk bcftools libopenblas-base wget nano

# Install R packages
RUN R -e "install.packages(c('stringr', 'dplyr', 'ggplot2', 'plotly', 'shinydashboard', 'shinyalert', 'DT', 'ape', 'htmlwidgets'))"

# Copy the app to the image
RUN rm -r /srv/shiny-server/*
COPY app.R /srv/shiny-server/
RUN mkdir /home/data
COPY mideca_streamdb.fasta /home/data/
COPY mifish_streamdb.fasta /home/data/

# Conda/Mamba installation
RUN cd tmp
RUN curl https://repo.anaconda.com/miniconda/Miniconda3-py310_23.3.1-0-Linux-x86_64.sh --output miniconda.sh
RUN bash miniconda.sh -bup /miniconda3
ENV PATH="/miniconda3/bin:$PATH"
RUN conda install -y -c conda-forge mamba

# Install base decona
RUN mkdir /home/github && \
    cd /home/github && \
    git clone https://github.com/ehill-iolani/decona.git
RUN sed -i -e "s/\r$//" /home/github/decona/install/install.sh
RUN bash /home/github/decona/install/install.sh
SHELL ["mamba", "run", "-n", "decona", "/bin/bash", "-c"]

# Add BLAST and medaka to decona
RUN mamba init && \
    mamba install -y -c bioconda blast=2.11.0 && \
    mamba install -y pandas=1.4.1 && \
    mamba install -y -c bioconda -c conda-forge bcftools=1.11 samtools=1.19.2 && \
    pip install medaka pyabpoa  && \
    echo "mamba activate decona" >> ~/.bashrc

# Clean up installation
RUN rm miniconda.sh

# Set working directory
WORKDIR /home/data

# Make /home/ writeable to all "users"
RUN chmod -R 777 /home/

# Make conda executable to all "users"
RUN chmod -R 777 /miniconda3/bin/*
