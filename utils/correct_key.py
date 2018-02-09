
import pandas as pd
import sys

def tidy(df):

    rowdicts = []
    for _, row in df.iterrows():
        row_dict = row.to_dict()

        baseQcovar = row_dict["baseQcovar"].split("|")

        for bqc in baseQcovar:
            row_dict_copy = row_dict.copy()
            row_dict_copy["baseQcovar"] = bqc

            rowdicts.append(row_dict_copy)

    ndf = pd.DataFrame(rowdicts)

    mdf = ndf.melt(value_vars="addBcovar addQcovar baseBcovar baseQcovar".split(),
                   id_vars="description N N_cases N_controls type version name".split())

    mdf = mdf.fillna("NA").drop_duplicates()

    mdf = mdf.rename(columns={"variable": "covar_type", "value": "covar_name"})

    return mdf["name	version	type	covar_type	covar_name	description	N	N_cases	N_controls".split()].reset_index(drop=True)


if __name__ == "__main__":

    key_file = sys.argv[1]

    df = pd.read_table(key_file, sep="\t")

    print(tidy(df).to_csv(sep="\t", index=False), file=sys.stdout, end="")
