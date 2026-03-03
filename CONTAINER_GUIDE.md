
# Genereate Container

### Create Dockerfile
  
`snakemake containerize > Dockerfile`

### Generate .def

manuell process. See [here](https://apptainer.org/docs/user/main/docker_and_oci.html#differences-and-limitations-vs-docker)

### Generate .sif

`apptainer build mycontainer.sif mycontainer.def`

# Run use-conda Workflow

### Open shell
`apptainer shell 
    --bind /etc/xdg/snakemake:/etc/xdg/snakemake 
    --bind /local/work
    --bind /groups/ds/resistance_cefiderocol/data/subsampled 
    --bind /local/work/julian/qc_pipeline/QC_pre_NextSeq:/workflow 
    Container_use_conda.sif`

### Activate snakemake env
`source /opt/conda/etc/profile.d/conda.sh`

`conda activate /conda-envs/30ef2e4d2f0920ae8c6c5118c93ed400`

### Change to workflow
`cd /workflow`

### Run Snakemake

`snakemake --cores all --use-conda`

# Run snakemake_env Workflow

### Open shell
`apptainer shell 
    --bind /etc/xdg/snakemake:/etc/xdg/snakemake 
    --bind /local/work
    --bind /groups/ds/resistance_cefiderocol/data/subsampled 
    --bind /local/work/julian/qc_pipeline/QC_pre_NextSeq:/workflow 
    Container_snakemake_env.sif`

### Activate snakemake env
`source /opt/conda/etc/profile.d/conda.sh`

`conda activate snakemake_env`

### Change to workflow
`cd /workflow`

### Run Snakemake

`snakemake --cores all`











apptainer shell 
    --overlay writable_overlay.img:rw
    --bind /etc/xdg/snakemake:/etc/xdg/snakemake 
    --bind /local/work
    --bind /groups/ds/resistance_cefiderocol/data/subsampled 
    --bind /local/work/julian/qc_pipeline/QC_pre_NextSeq:/workflow 
    Dockerfile_v5.sif


source /opt/conda/etc/profile.d/conda.sh

conda activate snakemake_env 


cd /workflow
snakemake     --cores all 

