.. HUNT-GWAS-Pipeline documentation master file, created by
   sphinx-quickstart on Fri Jan 26 13:20:19 2018.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to the HUNT-GWAS-Pipeline documentation!
================================================

.. image:: https://gitlab.com/huntgenes/gwas_pipeline/badges/master/build.svg

The HUNT-GWAS-Pipeline is a highly configurable pipeline that runs GWAS
analyses, plots the results and creates ample logs to help you understand which
steps were performed. No programming skill is required; the only thing the
pipeline needs is a filled out configuration file.

The HUNT-GWAS pipeline leverages the innovative Snakemake workflow language to
ensure reproducibility and scalability at no extra cost or additional complexity
for the user.

To ensure the reproducibility across machines and clusters, each step in the
workflow uses an accompanying Singularity-container with the necessary software
installed.

The workflow has ample tests to quickly allow for continual innovation and
contributions witout risking introducing new bugs.

.. toctree::
   :maxdepth: 2
   :caption: Getting started

   installation
   quick_start



Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
