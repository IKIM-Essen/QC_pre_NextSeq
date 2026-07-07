FROM condaforge/miniforge3:latest
LABEL io.github.snakemake.containerized="true"
LABEL io.github.snakemake.conda_env_hash="4f39cd98a09d47f2a68cabc77878eb9ec401443c18a8e58e86920d0410d8ddd5"
# Conda environment:
#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v7.6.0/bio/fastqc/environment.yaml
#   prefix: /conda-envs/136bb261dd776f5cb42eb4a819d921a5
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastqc =0.12.1
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/136bb261dd776f5cb42eb4a819d921a5
ADD https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v7.6.0/bio/fastqc/environment.yaml /conda-envs/136bb261dd776f5cb42eb4a819d921a5/environment.yaml
# Conda environment:
#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.4.2/bio/minimap2/aligner/environment.yaml
#   prefix: /conda-envs/8918f979401175af3caf79844d0ad783
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - minimap2 =2.30
#     - samtools =1.23.1
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/8918f979401175af3caf79844d0ad783
ADD https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.4.2/bio/minimap2/aligner/environment.yaml /conda-envs/8918f979401175af3caf79844d0ad783/environment.yaml
# Conda environment:
#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.4.2/bio/samtools/stats/environment.yaml
#   prefix: /conda-envs/f2d153d827f2b692d6500152b93cce9c
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - samtools =1.23.1
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/f2d153d827f2b692d6500152b93cce9c
ADD https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.4.2/bio/samtools/stats/environment.yaml /conda-envs/f2d153d827f2b692d6500152b93cce9c/environment.yaml
# Conda environment:
#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.8.0/bio/fastp/environment.yaml
#   prefix: /conda-envs/a44625d4cc5ac6471699c6dd12338e94
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - fastp =1.3.3
RUN mkdir -p /conda-envs/a44625d4cc5ac6471699c6dd12338e94
ADD https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.8.0/bio/fastp/environment.yaml /conda-envs/a44625d4cc5ac6471699c6dd12338e94/environment.yaml
# Conda environment:
#   source: https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.8.0/bio/multiqc/environment.yaml
#   prefix: /conda-envs/7e6193fb9857e98971b9985a859a93f4
#   channels:
#     - conda-forge
#     - bioconda
#     - nodefaults
#   dependencies:
#     - multiqc =1.34
#     - snakemake-wrapper-utils =0.8.0
RUN mkdir -p /conda-envs/7e6193fb9857e98971b9985a859a93f4
ADD https://raw.githubusercontent.com/snakemake/snakemake-wrappers/v9.8.0/bio/multiqc/environment.yaml /conda-envs/7e6193fb9857e98971b9985a859a93f4/environment.yaml
# Conda environment:
#   source: workflow/envs/kaiju_based.yaml
#   prefix: /conda-envs/4fc318cda14b65583913b27a213eaa1a
#   channels:
#     - conda-forge
#     - bioconda
#     - anaconda
#     - nodefaults
#   dependencies:
#     - kaiju
RUN mkdir -p /conda-envs/4fc318cda14b65583913b27a213eaa1a
COPY workflow/envs/kaiju_based.yaml /conda-envs/4fc318cda14b65583913b27a213eaa1a/environment.yaml
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

RUN conda env create --prefix /conda-envs/136bb261dd776f5cb42eb4a819d921a5 --file /conda-envs/136bb261dd776f5cb42eb4a819d921a5/environment.yaml && \
    conda env create --prefix /conda-envs/8918f979401175af3caf79844d0ad783 --file /conda-envs/8918f979401175af3caf79844d0ad783/environment.yaml && \
    conda env create --prefix /conda-envs/f2d153d827f2b692d6500152b93cce9c --file /conda-envs/f2d153d827f2b692d6500152b93cce9c/environment.yaml && \
    conda env create --prefix /conda-envs/a44625d4cc5ac6471699c6dd12338e94 --file /conda-envs/a44625d4cc5ac6471699c6dd12338e94/environment.yaml && \
    conda env create --prefix /conda-envs/7e6193fb9857e98971b9985a859a93f4 --file /conda-envs/7e6193fb9857e98971b9985a859a93f4/environment.yaml && \
    conda env create --prefix /conda-envs/4fc318cda14b65583913b27a213eaa1a --file /conda-envs/4fc318cda14b65583913b27a213eaa1a/environment.yaml && \
    conda env create --prefix /conda-envs/69ec81849c2d97f277bba0de2ad36894 --file /conda-envs/69ec81849c2d97f277bba0de2ad36894/environment.yaml && \
    conda env create --prefix /conda-envs/8bcfa5040c3f3fccb9990896a1047a49 --file /conda-envs/8bcfa5040c3f3fccb9990896a1047a49/environment.yaml && \
    conda env create --prefix /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400 --file /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400/environment.yaml && \
    conda clean --all -y
