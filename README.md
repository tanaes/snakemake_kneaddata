# snakemake_kneaddata
[Snakemake](http://snakemake.bitbucket.org) workflow for running 
[kneaddata](https://bitbucket.org/biobakery/kneaddata) for
basic quality control and human read removal from metagenomics data. 

## Background
This workflow is customized for execution on Barnacle, the Knight Lab's cluster
compute environment. It can be easily adapted to run in other environments by
modifying the cluster parameters in [cluster.json](./cluster.json) and
[launch.sh](./launch.sh.

## Installation
[install.sh](./install.sh) has a basic series of commands to set this workflow 
up in your Barnacle environment using Miniconda. This assumes you already have
Miniconda installed in your home directory; if not, you can install following
the directions [here](http://conda.pydata.org/miniconda.html).

To install this workflow, first make a new miniconda environment named
'kneaddata':

```
conda create -n kneaddata --yes python=3 pip pyyaml xlrd pandas jupyter notebook
```

Then enter the kneaddata environment and install kneaddata's dependencies:

```
source activate kneaddata
conda install -c bioconda --yes bowtie2 fastqc trimmomatic
```

Now install kneaddata, snakemake, and MultiQC through pip:

```
pip install kneaddata snakemake multiqc
```

Finally, we'll make a db directory in our conda environment for storing
kneaddata databases, and copy the demo db that comes with kneaddata to it:

```
mkdir -p $CONDA_ENV_PATH/share/kd_dbs
cp -r $CONDA_ENV_PATH/lib/python3.5/site-packages/kneaddata/tests/data/demo_bowtie2_db/ $CONDA_ENV_PATH/share/kd_dbs/demo
```

If you want to run against a full human genome database, run the following:

```
# kneaddata_database --download human bowtie2 $CONDA_ENV_PATH/share/kd_dbs/human
```

## Usage

This workflow requires the [Snakefile](./Snakefile) and a
[config.yaml](./config.yaml) that specifies filepaths and parameter values. 

An example file (and corresponding example data) has been provided in
[config_example_run.yaml](config_example_run.yaml). To execute the workflow
against these test data, make sure you are in the kneaddata environment,
navigate to the git repository directory, and execute:

```
snakemake --configfile config_example_run.yaml
```

To do a 'dry run' without executing, but simply list the commands that will be
run:

```
snakemake --configfile config_example_run.yaml -n
```

And to execute with the full cluster adaptation, run:

```
bash ./launch.sh ./ --configfile config_example_run.yaml
```

## Creating new config.yaml files

I have provided an iPython notebook [config_prep.ipynb](config_prep.ipynb) that
you can use toneasily generate config files for your own data. 

You will need the sequencing manifest provided by IGM, and the location of the
directory with the per-sample FASTQ files on Barnacle. The example ipynb can
then be modified to produce a config.yaml appropriate for your own samples. 