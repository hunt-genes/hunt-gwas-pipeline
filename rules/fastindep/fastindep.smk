
def get_king_kinship(config):

    if config.get("king_kinship", ""):
        king_kinship = config["king_kinship"]
    else:
        king_kinship = "{prefix}/data/king/kinship"

    return king_kinship


if remove_related:

    rule extract_unrelated:
        input:
            pheno = config["phenotype_file"],
            king_kinship = get_king_kinship(config)
            # fastindep_params = "{prefix}/data/extract_unrelated/params.txt"
        params:
            plink_binary_prefix = config["plink_binary_prefix"],
            prioritization = config["prioritize"],
            flip_priority = config["flip_priority"],
            unrelated_threshold = config["unrelated_threshold"],
            number_of_runs = config["number_of_runs"],
            seed = config["seed"]
        output:
            "{prefix}/data/extract_unrelated/{pheno}_unrelated.txt"
        conda:
            "../../envs/fastindep.yaml"
        benchmark:
            "{prefix}/benchmarks/extract_unrelated/{pheno}.benchmark"
        script:
            "../../scripts/01.5_ExtractUnrelated.R"
