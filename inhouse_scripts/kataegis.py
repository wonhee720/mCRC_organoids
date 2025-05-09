## 20210128 For NCH 

vcf_path = ""
df = scripts.tidydf(vcf_path).df

df.astype({"CHROM": 'str'})
df.astype({"POS": 'int'})
df[['context_3', 'substitution']] = df.query("VAR_TYPE == 'snp'").apply(lambda x: scripts.context_3(fasta, x['CHROM'], x['POS'], x['ALT']), axis=1, result_type='expand')
df.loc[df["VAR_TYPE"] != 'snp', 'substitution'] = 'others'
for chrom in chromosomes:
    df.loc[df['CHROM'] == chrom, 'logDIFF'] = np.log(df[df['CHROM'] == chrom]['POS'].diff())


fasta = pysam.FastaFile("/home/users/data/01_reference/human_g1k_v37/human_g1k_v37.fasta")


fig, ax = scripts.genome_figure()
fig.suptitle(f"{k} Kataegis Plot")
for i, chrom in enumerate(chromosomes[:-1]):
    df_tmp = df.query("CHROM == @chrom")[['CHROM', 'POS', 'substitution']]
    scripts.kataegis(df_tmp, ax[chrom], hue_column='substitution')
    if i != 0:
        try: 
            ax[chrom].get_legend().remove()
        except:
            pass