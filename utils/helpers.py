import re



def title(pk, w):

    trait_type = pk.loc[pk.name == w.pheno].type.iloc[0]
    trait = w.pheno

    if trait_type == "binary":

        n_cases = pk.loc[pk.name == w.pheno].N_cases.iloc[0]

        n_controls = pk.loc[pk.name == w.pheno].N_controls.iloc[0]

        title = "".join([str(s) for s in ["Phenotype: ", trait, " Type: ",
                                        trait_type + "\n(", n_cases, " cases versus ",
                                        n_controls, " controls)"]])

    elif trait_type == "quantitative":

        n = pk.loc[pk.name == w.pheno].N.iloc[0]

        title = "".join([str(s) for s in ["Phenotype: ", trait, " Type: ",
                                          trait_type + ", N=", n]])


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
