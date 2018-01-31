from subprocess import run

import pandas as pd


def run_dag(targets, configfile, snakefile="Snakefile", config="", extras="", dryrun=True): #, create_env=False

    # abs_configfile = abspath(configfile)
    # abs_sample_sheet = abspath(sample_sheet)

    # print(abs_configfile)
    # print(abs_sample_sheet)

    # with TemporaryDirectory() as tempdir:
    config = "--config " + config

    dryrun = "n" if dryrun else ""

    cmd = ("snakemake -s {snakefile} --configfile {configfile}"
            " -{dryrun}p {targets} -F {extras} {config}").format(**locals())

    print(cmd)

    return run(cmd, shell=True).returncode

def fetch_main_targets(snakefile="Snakefile"):
    "Get all targets in main Snakefile"

    r = re.compile("^rule (.*):")

    targets = []
    for line in open(snakefile):
        match = r.search(line)
        if match:
            targets.append(match.group(1))

    return targets
