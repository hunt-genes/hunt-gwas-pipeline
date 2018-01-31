
Installation
================================================

The HUNT-GWAS Pipeline requires Snakemake and the conda package manager.
Singularity is optional, but recommended. This page will show you how to install
them.

First we will install the conda pacakge manager. Go to
https://www.anaconda.com/download/ and select the Python 3+ version.

Then follow the installation instructions.

Next, we will need to install snakemake. We will do this using anaconda.

.. code-block:: bash

   # install snakemake from the channel bioconda
   conda install -c bioconda snakemake

Finally you can install the HUNT-GWAS Pipeline with the command

.. code-block:: bash

   git clone git@gitlab.com:huntgenes/gwas_pipeline.git

Singularity
~~~~~~~~~~~

Singularity is a way to package software and the operating system it is run on
to ensure complete reproducibility. Using it with the HUNT-GWAS pipeline is
highly recommended. The installation-instructions are here:
http://singularity.lbl.gov/install-linux
