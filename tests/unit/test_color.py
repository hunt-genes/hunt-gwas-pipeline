import pytest

import pandas as pd
from io import StringIO

@pytest.fixture
def df_gwas_catalogue():

    c = """DATE ADDED TO CATALOG	PUBMEDID	FIRST AUTHOR	DATE	JOURNAL	LINK	STUDY	DISEASE/TRAIT	INITIAL SAMPLE SIZE	REPLICATION SAMPLE SIZE	REGION	CHR_ID	CHR_POS	REPORTED GENE(S)	MAPPED_GENE	UPSTREAM_GENE_ID	DOWNSTREAM_GENE_ID	SNP_GENE_IDS	UPSTREAM_GENE_DISTANCE	DOWNSTREAM_GENE_DISTANCE	STRONGEST SNP-RISK ALLELE	SNPS	MERGED	SNP_ID_CURRENT	CONTEXT	INTERGENIC	RISK ALLELE FREQUENCY	P-VALUE	PVALUE_MLOG	P-VALUE (TEXT)	OR or BETA	95% CI (TEXT)	PLATFORM [SNPS PASSING QC]	CNV
2009-09-28	18403759	Ober C	2008-04-09	N Engl J Med	www.ncbi.nlm.nih.gov/pubmed/18403759	Effect of variation in CHI3L1 on serum YKL-40 level, risk of asthma, and lung function.	YKL-40 levels	632 Hutterite individuals	443 European ancestry cases, 491 European ancestry controls, 206 European ancestry individuals	1q32.1	1	203186754	CHI3L1	CHI3L1			1116			rs4950928-G	rs4950928	0	4950928	upstream_gene_variant	0	0.29	1E-13	13.0		0.3	[NR] ng/ml decrease	Affymetrix [290325]	N
2008-06-16	18369459	Liu Y	2008-04-04	PLoS Genet	www.ncbi.nlm.nih.gov/pubmed/18369459	A genome-wide association study of psoriasis and psoriatic arthritis identifies new disease loci.	Psoriasis	218 European ancestry cases, 519 European ancestry controls	1,153 European ancestry cases, 1,217 European ancestry controls	13q14.11	13	39776775	COG6	COG657511			rs7993214-?	rs7993214	0	7993214	intron_variant	0	0.65	2E-6	5.698970004336019		1.41	[1.22-1.61]	Illumina [305983]	N
2008-06-16	18385676	Amos CI	2008-04-03	Nat Genet	www.ncbi.nlm.nih.gov/pubmed/18385676	Genome-wide association scan of tag SNPs identifies a susceptibility locus for lung cancer at 15q25.1.	Lung cancer	1,154 European ancestry cases, 1,137 European ancestry controls	2,724 European ancestry cases, 3,694 European ancestry controls	15q25.1	15	78513681	CHRNA5, CHRNA3, PSMA4, LOC123688	HYKK			123688			rs8034191-G	rs8034191	0	8034191	intron_variant	0	NR	3E-18	17.522878745280337		1.3	[1.15-1.47]	Illumina [317498]	N
2008-06-16	18385676	Amos CI	2008-04-03	Nat Genet	www.ncbi.nlm.nih.gov/pubmed/18385676	Genome-wide association scan of tag SNPs identifies a susceptibility locus for lung cancer at 15q25.1.	Lung cancer	1,154 European ancestry cases, 1,137 European ancestry controls	2,724 European ancestry cases, 3,694 European ancestry controls	1q23.2	1	159711078	CRP	CRPP1 - CRP	171422	1401		5482	1211	rs2808630-G	rs2808630	0	2808630	downstream_gene_variant	1	NR	7E-6	5.154901959985743		1.22	[1.10-1.35]	Illumina [317498]	N
2008-06-16	18385676	Amos CI	2008-04-03	Nat Genet	www.ncbi.nlm.nih.gov/pubmed/18385676	Genome-wide association scan of tag SNPs identifies a susceptibility locus for lung cancer at 15q25.1.	Lung cancer	1,154 European ancestry cases, 1,137 European ancestry controls	2,724 European ancestry cases, 3,694 European ancestry controls	3q28	3	190632672	IL1RAP	IL1RAP			3556			rs7626795-G	rs7626795	0	7626795	intron_variant	0	NR	8E-6	5.096910013008056		1.16	[1.05-1.28]	Illumina [317498]	N
2008-06-16	18385738	Hung RJ	2008-04-03	Nature	www.ncbi.nlm.nih.gov/pubmed/18385738	A susceptibility locus for lung cancer maps to nicotinic acetylcholine receptor subunit genes on 15q25.	Lung cancer	1,926 European ance other ancestry cases, 2,522 European and other ancestry controls	332 European ancestry cases, 462 European ancestry controls, 2,181 European and other ancestry cases, 4,290 European and other ancestry controls	15q25.1	15	78513681	CHRNA5, CHRNA3, PSMA4, LOC123688, CHRNB4	HYKK			123688			rs8034191-C	rs8034191	0	8034191	intron_variant	0	0.34	5E-20	19.30102999566398		1.3	[1.23-1.37]	Illumina [310023]	N
2008-09-16	18385739	Thorgeirsson TE	2008-04-03	Nature	www.ncbi.nlm.nih.gov/pubmed/18385739	A variant associated with nicotine dependence, lung cancer and peripheral arterial disease.Nicotine dependence	10,995 European ancestry individuals	4,848 European ancestry individuals	15q25.1	15	78601997	CHRNA5, CHRNA3, CHRNB4	CHRNA3			1136			rs1051730-T	rs1051730	0	1051730	synonymous_variant	0	0.35	6E-20	19.221848749616356		0.1	[0.08-0.12] cigarettes per day increase	Illumina [306207]	N
2008-06-16	18372901	Tenesa A	2008-03-30	Nat Genet	www.ncbi.nlm.nih.gov/pubmed/18372901	Genome-wide association scan identifies a colorectal cancer susceptibility locus on 11q23 and replicates risk loci at 8q24 and 18q21.	Colorectal cancer	981 European ancestry cases, 1,002 European ancestry controls	10,287 European ancestry cases, 10,401 European ancestry controls, 4,400 Japanese ancestry cases, 3,179 Japanese ancestry controls, 1,789 cases, 1,771 controls	8q24.21	8	127412547	POU5FIP1, HsG57825, DQ515897	CASC8			727677		rs7014346-A	rs7014346	0	7014346	intron_variant	0	0.18	9E-26	25.045757490560675		1.19	[1.15-1.23]	Illumina [541628]	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	8q13.3	8	71601993	MSC	LOC105375892 - LOC107986890	105375892	107986890		46818	129752	rs12679556-T	rs12679556	0	12679556	intergenic_variant	1	0.7074	5E-8	7.301029995663981				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	12p11.23	12	26318431	ITPR2, SSPN	LOC105369705 - LOC107984482	105369705	107984482		7967	4489	rs10842707-T	rs10842707	0	10842707	intron_variant	1	0.2246	5E-6	5.301029995663981				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	12q13.13	12	53948900	HOXC13	HOXC13 - LOC105369775	3229	105369775	2356	4485	rs1443512-A	rs1443512	0	1443512	downstream_gene_variant	1	0.2662	5E-7	6.301029995663981				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	12q13.13	12	53955989	HOXC12	HOXC12, LOC105369775			3228, 105369775			rs10783615-A	rs10783615	0	10783615	intron_variant	0	0.7973	2E-10	9.698970004336019				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	12q24.31	12	123955563	CCDC92	CCDC92			80212			rs4765219-A	rs4765219	0	4765219	intron_variant	0	0.3303	2E-7	6.698970004336019				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	15q13.3	15	31416060	KLF13	KLF13			51621			rs8042543-T	rs8042543	0	8042543	intron_variant	0	0.2071	7E-6	5.154901959985743				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-08-10	28443625	Justice AE	2017-04-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28443625	Genome-wide meta-analysis of 241,258 adults accounting for smoking behaviour identifies novel loci for obesity traits.	Waist-to-hip ratio adjusted for BMI (joint analysis main effects and smoking interaction)	95,362 European ancestry women, 62,085 European ancestry men, 5,829 European ancestry individuals, 10,500 African American/Afro-Caribbean ancestry women, 2,706 African American/Afro-Caribbean ancestry men, 1,030 Indian Asian ancestry women, 7,648 Indian Asian ancestry men, 1,793 Filipino ancestry women, 2,944 Hispanic/Latino ancestry women, 1,764 Hispanic/Latino ancestry men	21,496 European ancestry women, 24,385 European ancestry men, 118,364 European ancestry individuals, 2,326 African American/Afro-Caribbean ancestry women, 855 African American/Afro-Caribbean ancestry men	20q13.12	20	46930192	EYA2	EYA2			2139			rs6090583-A	rs6090583	0	6090583	intron_variant	0	0.495	3E-7	6.522878745280337				Affymetrix, Illumina, Perlegen [up to 2800000] (imputed)	N
2017-09-20	28548082	Southam L	2017-05-26	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28548082	Whole genome sequencing and imputation in isolated populations identify genetic associations with medically-relevant complex traits.	Fasting blood glucose adjusted for BMI	up to 1,476 Mylopotamos (founder/genetic isolate) individuals, up to 1,737 Pomak (founder/genetic isolate) individuals	NA	20p12.2	20	10453882	SLX4IP	SLX4IP			128710			rs6131100-A	rs6131100	0	6131100	intron_variant	0	0.037	1E-8	8.0(Pomak)	0.79	[0.52-1.06] unit decrease	Illumina [at least 13541454] (imputed)	N
2017-10-06	28598419	Hong X	2017-06-09	Nat Commun	www.ncbi.nlm.nih.gov/pubmed/28598419	Genome-wide approach identifies a novel gene-maternal pre-pregnancy BMI interaction on preterm birth.	Preterm birth (maternal effect) (maternal pre-pregnancy BMI interaction)	329 African American normal weight preterm birth mothers, 514 African American normal weight term birth mothers, 463 African American overweight or obese preterm birth mothers, 580 African American overweight or obese term birth mothers	NA	1p22.3	1	86022231	COL24A1	COL24A1			255631			rs11161721-A	rs11161721	0	11161721	intron_variant	0	.191	4E-9	8.397940008672037				Affymetrix, Illumina [9929081] (imputed)	N"""

    return pd.read_table(StringIO(c), low_memory=False)


@pytest.fixture
def parsed_df():

    c = """CHR_ID CHR_POS TRAIT
1 750 BMI
1 1000 HEIGHT"""

@pytest.fixture
def expected_result():

    c = """CHROM;START;END;COL;LEGENDTEXT
1;500;1000;BLUE;'BMI SNPS'"""

    return pd.read_table(StringIO(c), low_memory=False)


def manhattan_plot_regions(df, slack=500000):

    new_dfs = []
    for trait, trait_df in df.groupby("TRAIT"):
        tdf = trait_df

        start = tdf.CHR_POS - slack
        end = tdf.CHR_POS + slack



        new_dfs.append(trait_df)

    return pd.concat(new_dfs)










@pytest.mark.unit
def test_gwas_catalogue(parsed_df):

    pass
    # result = manhattan_plot_regions(parsed_df)

    # assert expected_result.equals(parsed_df)
