# how to install this workflow

# make a conda env for kneaddata
conda create -n kneaddata python=3 pip snakemake pyyaml xlrd pandas jupyter notebook

# install kneaddata and dependencies
source activate kneaddata

conda install -c bioconda bowtie2 fastqc

pip install kneaddata

# download the kneaddata human genome database
