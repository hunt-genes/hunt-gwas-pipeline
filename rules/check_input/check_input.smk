# rule check_input:
#     input:
#         pheno = config["phenotype_file"],
#         key = config["key_file"]
#     output:
#         expand("{prefix}/data/check_input/{pheno}_pheno.txt", pheno=phenotypes, prefix=prefix),
#         expand("{prefix}/data/check_input/{pheno}_summary.txt", pheno=phenotypes, prefix=prefix),
#         checked = "{prefix}/data/check_input/checked.txt".format(prefix=prefix)
#     params:
#         outdir = "{prefix}/data/check_input/".format(prefix=prefix),
#         normalize = True
#     conda:
#         f"{base_prefix}/envs/check_input.yaml"
#     benchmark:
#         f"{prefix}/benchmarks/check_input/benchmark"
#     wrapper:
#         f"file://{base_prefix}/wrappers/check_input_py"


rule check_input_py:
    input:
        pheno = config["phenotype_file"],
        key = config["key_file"]
    output:
        f"{prefix}/data/check_input/{{pheno}}_pheno.txt",
        f"{prefix}/data/check_input/{{pheno}}_summary.txt",
        checked = "{prefix}/data/check_input/{{pheno}}_checked.txt".format(prefix=prefix)
    params:
        outdir = "{prefix}/data/check_input/".format(prefix=prefix),
        normalize_quantitative_traits = True
    conda:
        f"{base_prefix}/envs/check_input.yaml"
    benchmark:
        f"{prefix}/benchmarks/check_input/benchmark"
    wrapper:
        f"file://{base_prefix}/wrappers/check_input_py"
