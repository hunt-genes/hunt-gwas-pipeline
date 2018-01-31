# import sys

# from scipy.stats import chi2
# from scipy.special import chdtri
# import pandas as pd
# import numpy as np

# imputed_infile = snakemake.input.imputed
# # summary_infile = snakemake.input.summary

# df = pd.read_table(imputed_infile, header=0)
# # pheno = pd.read_table(summary_infile, header=0)

# trait = snakemake.wildcards["pheno"]

# binary = snakemake.params.trait_type

# analysis = snakemake.params.analysis

# maf_threshold = snakemake.params.maf_threshold

# # names of columns in SAIGE analysis #
# pcol = "p.value"
# ycol = "log10P_GC"
# chicol = "CHISQ"
# betacol = "BETA"
# secol = "SE"
# chrcol = "CHR"
# poscol = "POS"
# idcol = "ID"
# refcol = ""
# altcol = ""
# accol = ""
# maccol = ""
# afcol = "AF"
# af0col = ""
# af1col = ""

# # remove rows with no SE
# df = df.loc[df[secol] > 0]

# new_ps = -np.log10(df[pcol])

# no_p_value = df[pcol] == 0

# if len(no_p_value) > 0 and chicol in df:
#     new_ps.loc[no_p_value] = np.log(1 - chi2.cdf(new_ps.loc[no_p_value], df=1))
# elif len(no_p_value) > 0 and chicol in df:
#     new_ps.loc[no_p_value] = pd.Series([sys.float_info.min] * len(no_p_value))

# maf = df[afcol]
# maf.loc[maf > 0.5] = 1 - maf[maf > 0.5]

# chisq = chdtri(1, new_ps)
# print(chisq)

# min_maf = maf.min()


# def log10toP(log10_p_values):

#     abs_log10_p_values = abs(log10_p_values)

#     assert abs_log10_p_values.isnull().sum() == 0

#     if abs_log10_p_values > 300:

#         part1 = abs_log10_p_values // (100 * 100)
#         part2 = abs_log10_p_values - part1

#         if part2 != 0:
