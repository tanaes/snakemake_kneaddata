# how to install this workflow

# make a conda env for kneaddata
conda create -n kneaddata --yes python=3 pip snakemake pyyaml xlrd pandas jupyter notebook

# install kneaddata and dependencies
source activate kneaddata

conda install -c bioconda --yes bowtie2 fastqc trimmomatic

pip install kneaddata

# download the kneaddata human genome database
mkdir -p $CONDA_ENV_PATH/share/kd_dbs
cp -r $CONDA_ENV_PATH/lib/python3.5/site-packages/kneaddata/tests/data/demo_bowtie2_db/ $CONDA_ENV_PATH/share/kd_dbs/demo
# kneaddata_database --download human bowtie2 $CONDA_ENV_PATH/share/kd_dbs/human