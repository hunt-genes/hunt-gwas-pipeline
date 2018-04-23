
def maf_threshold(w):

    binary = pk.loc[pk.name == w.pheno].type.iloc[0] == "binary"

    min_mac = config["min_minor_allele_count"]

    if binary:
        n_controls = pk[pk.name == w.pheno].N_controls.iloc[0]
        n_cases = pk.loc[pk.name == w.pheno].N_cases.iloc[0]
        maf = min_mac / (2 * min(n_cases, n_controls))
    else:
        n = pk.loc[pk.name == w.pheno].N.iloc[0]
        maf = min_mac / (2 * n)

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


# rule prepare_for_plotting_saige:
#     input:
#         "{prefix}/data/collect_imputed_files/{pheno}.imputed.results"
#     output:
#         "{prefix}/data/imputed/{pheno}.final_results.txt"
#     run:
        # columns = "SNP	CHR     CM      POS     EFFECT_ALLELE   ALT_ALLELE      AC      AF      N       BETA    SE      Tstat   p.value p.value.NA      Is.SPA.converge varT    varTstar".split()

        # usecols = "CHR SNP POS p.value AF BETA SE".split()
        # df = pd.read_table(input[0], header=0)

        # chisq = (df.BETA / df.SE) ** 2
