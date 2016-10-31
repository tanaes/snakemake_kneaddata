# how to install this workflow

# make a conda env for kneaddata
conda create -n kneaddata --yes python=3 pip pyyaml xlrd pandas jupyter notebook

# install kneaddata and dependencies
source activate kneaddata

conda install -c bioconda --yes bowtie2 fastqc trimmomatic

pip install kneaddata snakemake multiqc

# download the kneaddata human genome database
mkdir -p $CONDA_ENV_PATH/share/kd_dbs
cp -r $CONDA_ENV_PATH/lib/python3.5/site-packages/kneaddata/tests/data/demo_bowtie2_db/ $CONDA_ENV_PATH/share/kd_dbs/demo
# kneaddata_database --download human bowtie2 $CONDA_ENV_PATH/share/kd_dbs/human


# Download metaphlan2
wget https://bitbucket.org/biobakery/metaphlan2/get/default.zip
unzip default.zip 'biobakery-metaphlan2*/db_v20/*' 'biobakery-metaphlan2*/metaphlan2.py'
mv biobakery-metaphlan2*/metaphlan2.py $CONDA_ENV_PATH/bin/.
mv biobakery-metaphlan2*/db_v20 $CONDA_ENV_PATH/bin/metaphlan_db