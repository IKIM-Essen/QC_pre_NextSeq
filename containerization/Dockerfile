FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="bc8c7110aaf156d82bd725d0e6ccaa1d4fd89ce15a97293689faa4eb5823cfa7"

# Step 2: Retrieve conda environments

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/fastp/environment.yaml
#   prefix: /conda-envs/51cfe917ef083f45a131d4bed3b30579
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastp =0.23.4
RUN mkdir -p /conda-envs/51cfe917ef083f45a131d4bed3b30579
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/fastp/environment.yaml /conda-envs/51cfe917ef083f45a131d4bed3b30579/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/fastqc/environment.yaml
#   prefix: /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/24b8923f8e4abe077ffe95b01bfc1652
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/fastqc/environment.yaml /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/minimap2/aligner/environment.yaml
#   prefix: /conda-envs/b6e7a40ecb21f1430a0e790a9b12d13e
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - minimap2 =2.26
#     - samtools =1.19
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/b6e7a40ecb21f1430a0e790a9b12d13e
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/minimap2/aligner/environment.yaml /conda-envs/b6e7a40ecb21f1430a0e790a9b12d13e/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/multiqc/environment.yaml
#   prefix: /conda-envs/14870f1e101622320c0226eacb79a0ff
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.19
RUN mkdir -p /conda-envs/14870f1e101622320c0226eacb79a0ff
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/multiqc/environment.yaml /conda-envs/14870f1e101622320c0226eacb79a0ff/environment.yaml

# Conda environment:
#   source: https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/samtools/stats/environment.yaml
#   prefix: /conda-envs/02e79c5681016224e152a62d224d1b94
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.19
#     - snakemake-wrapper-utils =0.6.2
RUN mkdir -p /conda-envs/02e79c5681016224e152a62d224d1b94
ADD https://github.com/snakemake/snakemake-wrappers/raw/v3.3.3/bio/samtools/stats/environment.yaml /conda-envs/02e79c5681016224e152a62d224d1b94/environment.yaml

# Conda environment:
#   source: workflow/envs/kraken_based.yaml
#   prefix: /conda-envs/f338784d62bf13686a59e1955c303938
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - nodefaults
#   dependencies:
#     - kraken2 = 2.1.6
#     - krakentools = 1.2.1
#     - bracken=3.1
RUN mkdir -p /conda-envs/f338784d62bf13686a59e1955c303938
COPY workflow/envs/kraken_based.yaml /conda-envs/f338784d62bf13686a59e1955c303938/environment.yaml

# Conda environment:
#   source: workflow/envs/python.yaml
#   prefix: /conda-envs/69ec81849c2d97f277bba0de2ad36894
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#   dependencies:
#     - pandas = 2.1.4
#     - altair = 5.2.0
#     - numpy = 1.26.3
#     - distinctipy = 1.3.4
#     - matplotlib = 3.8.2
RUN mkdir -p /conda-envs/69ec81849c2d97f277bba0de2ad36894
COPY workflow/envs/python.yaml /conda-envs/69ec81849c2d97f277bba0de2ad36894/environment.yaml

# Conda environment:
#   source: workflow/envs/rbt.yaml
#   prefix: /conda-envs/8bcfa5040c3f3fccb9990896a1047a49
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - rust-bio-tools = 0.42.2
RUN mkdir -p /conda-envs/8bcfa5040c3f3fccb9990896a1047a49
COPY workflow/envs/rbt.yaml /conda-envs/8bcfa5040c3f3fccb9990896a1047a49/environment.yaml

# Conda environment:
#   source: workflow/envs/snakemake.yaml
#   prefix: /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - snakemake>=9.1
#     - pip:
#       - snakemake-storage-plugin-fs
RUN mkdir -p /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400
COPY workflow/envs/snakemake.yaml /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400/environment.yaml

# Conda environment:
#   source: workflow/envs/unix.yaml
#   prefix: /conda-envs/2e0b9e54d7d4d87567c5ffb2cc26709c
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - wget
RUN mkdir -p /conda-envs/2e0b9e54d7d4d87567c5ffb2cc26709c
COPY workflow/envs/unix.yaml /conda-envs/2e0b9e54d7d4d87567c5ffb2cc26709c/environment.yaml

# Step 3: Generate conda environments

RUN conda env create --prefix /conda-envs/51cfe917ef083f45a131d4bed3b30579 --file /conda-envs/51cfe917ef083f45a131d4bed3b30579/environment.yaml && \
    conda env create --prefix /conda-envs/24b8923f8e4abe077ffe95b01bfc1652 --file /conda-envs/24b8923f8e4abe077ffe95b01bfc1652/environment.yaml && \
    conda env create --prefix /conda-envs/b6e7a40ecb21f1430a0e790a9b12d13e --file /conda-envs/b6e7a40ecb21f1430a0e790a9b12d13e/environment.yaml && \
    conda env create --prefix /conda-envs/14870f1e101622320c0226eacb79a0ff --file /conda-envs/14870f1e101622320c0226eacb79a0ff/environment.yaml && \
    conda env create --prefix /conda-envs/02e79c5681016224e152a62d224d1b94 --file /conda-envs/02e79c5681016224e152a62d224d1b94/environment.yaml && \
    conda env create --prefix /conda-envs/f338784d62bf13686a59e1955c303938 --file /conda-envs/f338784d62bf13686a59e1955c303938/environment.yaml && \
    conda env create --prefix /conda-envs/69ec81849c2d97f277bba0de2ad36894 --file /conda-envs/69ec81849c2d97f277bba0de2ad36894/environment.yaml && \
    conda env create --prefix /conda-envs/8bcfa5040c3f3fccb9990896a1047a49 --file /conda-envs/8bcfa5040c3f3fccb9990896a1047a49/environment.yaml && \
    conda env create --prefix /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400 --file /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400/environment.yaml && \
    conda env create --prefix /conda-envs/2e0b9e54d7d4d87567c5ffb2cc26709c --file /conda-envs/2e0b9e54d7d4d87567c5ffb2cc26709c/environment.yaml && \
    conda clean --all -y
