import pandas as pd

if not config:
    configfile: "config.yaml"

dosage_column_names = "dosage_chr_col dosage_snp_id_col dosage_pos_col dosage_ref_allele_col dosage_alt_allele_col"
for k in dosage_column_names.split():
    assert config[k] in config["dosage_file_colnames_skip"], "Mismatch in dosage_file_colnames_skip and dosage_column names in config"


# pheno keys
pk = pd.read_table(config["key_file"])

phenotypes = pk.name
pheno_trait_type = {k: v for k, v in zip(phenotypes, pk["type"])}

dosage_df = pd.read_table(config["sample_sheet"])
dosages = {k: v for k, v in zip(dosage_df.Name, dosage_df.Path)}

prefix = config["intermediate_prefix"]
# base_prefix = "/home/endrebak/code/gwas"

remove_related = config["remove_related"]

wildcard_constraints:
    pheno = "|".join(phenotypes),
    chunk = "|".join(dosages.keys())

id_column = config["id_column"]
num_markers = config["number_of_markers"]

pheno_all = list(phenotypes)

# "check_input/check_input"
to_include = ["fastindep/fastindep", "saige/step1",
              "saige/step2", "collect_chunks/collect_chunks",
              "prepare_for_plotting/prepare_for_plotting", "qqplot/qqplot",
              "manhattan/manhattan"]


for rule in to_include:
    include: "rules/{rule}.smk".format(rule=rule)


rule all:
    input:
        expand("{prefix}/data/manhattan/{pheno}.imputed_Manhattan.png", prefix=prefix, pheno=pheno_all),
        # expand("{prefix}/data/qqplot/{pheno}.QQ_Plot.png", prefix=prefix, pheno=pheno_all),
        # expand("{prefix}/data/manhattan/{pheno}.imputed_Manhattan.png", prefix=prefix, pheno=pheno_all)


rule qq:
    input:
        expand("{prefix}/data/qqplot/{pheno}.QQ_Plot.png", prefix=prefix, pheno=pheno_all)


rule manhattan:
    input:
        expand("{prefix}/data/manhattan/{pheno}.imputed_Manhattan.png", prefix=prefix, pheno=pheno_all)


rule imputation:
    input:
        expand("{prefix}/data/collect_imputed_files/{pheno}.imputed.results", prefix=prefix, pheno=pheno_all)