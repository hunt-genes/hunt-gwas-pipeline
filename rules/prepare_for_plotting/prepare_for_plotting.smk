
def maf_threshold(w):

    binary = pk.loc[pk.name == w.pheno].type.iloc[0] == "binary"
    n_cases = pk.loc[pk.name == w.pheno].N_cases.iloc[0]

    min_mac = config["min_minor_allele_count"]

    if binary:
        n_controls = pk[pk.name == w.pheno].N_controls.iloc[0] if binary else 0
        maf = min_mac / (2 * min(n_cases, n_controls))
    else:
        maf = min_mac / (2 * n_cases)

    return maf


# Cannot be easily multithreaded due to candidate regions
# at least need to split into chromos first, but that might be as time-consuming
rule prepare_for_plotting:
    input:
        imputed = "{prefix}/data/collect_imputed_files/{pheno}.imputed.results"
    output:
        final_results = "{prefix}/data/imputed/{pheno}.final_results.txt",
        candidate_regions = "{prefix}/data/imputed/{pheno}.imputed.candidate_regions"
    params:
        analysis = "SAIGE", # only one implemented yet
        maf_threshold = maf_threshold,
        # trait_type = lambda w: "binary" if (pk.loc[pk.name == w.pheno].type.iloc[0] == "binary") else "quantitative"
    benchmark:
        "{prefix}/benchmarks/prepare_for_plotting/{pheno}.benchmark"
    script:
        "./../../scripts/04_prepareResults.R"
