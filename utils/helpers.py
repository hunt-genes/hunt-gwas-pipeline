import re


def trait_type(pk, w):

    binary = pk[w.pheno].type == "binary"
    trait_type = "binary" if binary else "quantitative"

    return trait_type


def title(pk, w):

    binary = pk.loc[pk.name == w.pheno].type.iloc[0] == "binary"
    trait_type = "binary" if binary else "quantitative"
    n_cases = pk.loc[pk.name == w.pheno].N_cases.iloc[0]
    trait = w.pheno

    n_controls = pk.loc[pk.name == w.pheno].N_controls.iloc[0] if binary else 0

    title = "".join([str(s) for s in ["Phenotype: ", trait, " Type: ",
                                      trait_type + "\n(", n_cases, " cases versus ",
                                      n_controls, " controls)"]])

    return title


def fetch_main_targets(snakefile="Snakefile"):
    "Get all targets in main Snakefile"

    r = re.compile("^rule (.*):")

    targets = []
    for line in open(snakefile):
        match = r.search(line)
        if match:
            targets.append(match.group(1))

    return targets
