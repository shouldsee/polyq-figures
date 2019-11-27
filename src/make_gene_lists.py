import pymisca.ext as pyext
import src.rnaseq_figure
import src.util as _util
df = pyext.OrderedDict()
df["signature_targets"] = (src.rnaseq_figure.job['signature_targets'])
df["chipseq_targets_genes"] = pyext.readData("OUTPUT/chipseq_targets_genes_job.peak_list.csv")['feat_acc'].unique()
# df = pyext.pd.DataFrame(map(pyext.pd.Series(df.items())))
df = pyext.pd.DataFrame.from_dict(df,orient='index').T
df.to_csv(_util._get_output_file("OUTPUT/gene_lists_dataframe.csv"),index=0)