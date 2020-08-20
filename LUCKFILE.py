# %%python2
# from rnaseq_figure import job as template
from luck.shorts import RNS, TSSR, NCR, LSC
RULE = TSSR
ns = RNS()

from src.util import _get_output_file, Inputer
import os
SRCDIR = os.path.realpath(os.getcwd())


RULE.M(ns, "src/util.py")
if 1:
    _ = '''
    #####################################
    #### Making RNASEQ boxplots  ########
    '''
    inputs = Inputer()
    input_rnaseq = inputs._get_file('/home/feng/envs/Fig_POLYQ/rnaseq.pk')
    input_datasets_meta = inputs._get_file('/home/feng/static/results/0318-makeRNA-polyQ/mcurr0.csv')
    input_markers_df = inputs._get_file(f'{SRCDIR}/src/key_ath.csv',raw=1)
    for FNAME in inputs:
        RULE.MWF(ns, FNAME, None)

    SCRIPT = "src/rnaseq_figure_2.py"
    RULE.MWF(ns, SCRIPT, None)

    OUTS = []
    for FNAME in [
        'src/BOXPLOT_1_OX-lines.tsv',
        'src/BOXPLOT2_mutants.tsv',
        'src/BOXPLOT3_Qlines.tsv',
    ]:
        # job = template.copy()
        FNAME = FNAME 
        OFNAME = "OUTPUT/" + FNAME.replace('/','_')+'.svg'
        OFNAME = _get_output_file(OFNAME)
        RULE.MWF(ns, FNAME, None)
        RULE.MWF(ns, OFNAME, ' '.join([FNAME] + [SCRIPT, input_rnaseq, input_datasets_meta, input_markers_df]),  f"\
            \npython2 {SCRIPT} \
            --input_rnaseq {input_rnaseq} \
            --input_datasets_meta {input_datasets_meta} \
            --input_markers_df {input_markers_df} \
            --FNAME {FNAME} \
            --OFNAME {OFNAME} \
            --TARGET visualise_boxplot \
            --run \
            ")
        OUTS.append(OFNAME)
    NCR.M(ns, "boxplots", ' '.join(OUTS))

    _ = '''
    #####################################
    #### Making RNASEQ heatmaps  ########
    '''
    OUTS = []
    for FNAME in [
        'src/FIG-S6_HEATMAP1_Q-lines.tsv',
        'src/FIG-S5_HEATMAP2_genotypes.tsv',
        'src/FIG-S7_HEATMAP3.tsv',
    ]:
        # job = template.copy()
        FNAME = FNAME 
        OFNAME = "OUTPUT/" + FNAME.replace('/','_')+'.png'
        OFNAME = _get_output_file(OFNAME)

        RULE.MWF(ns, FNAME, None)
        RULE.MWF(ns, OFNAME, 
            ' '.join([FNAME] + [SCRIPT, input_rnaseq, input_datasets_meta, input_markers_df]),  f"\
            \npython2 {SCRIPT} \
            --input_rnaseq {input_rnaseq} \
            --input_datasets_meta {input_datasets_meta} \
            --input_markers_df {input_markers_df} \
            --FNAME {FNAME} \
            --OFNAME {OFNAME} \
            --TARGET visualise_heatmap \
            --run \
            ")
        OUTS.append(OFNAME)
    NCR.M(ns, "heatmaps", ' '.join(OUTS))


if 1:
    _ = '''
    #############################
    ### Make chipseq targets ####
    '''
    SCRIPT = "src/make_chipseq_targets.py"
    OUTPUTS= LSC(f'python2 {SCRIPT} --print-outputs').decode().replace("\n"," ")
    INPUTS =  LSC(f'python2 {SCRIPT} --print-inputs').decode().replace("\n",' ')
    INPUTS += " "+ SCRIPT
    # from pprint import pprint
    # pprint([INPUTS,OUTPUTS])
    # import pdb; pdb.set_trace()
    [ RULE.MWF(ns, INPUT, None) for INPUT in INPUTS.split()]
    # assert 0,OUTPUTS
    RULE.MWF(ns, OUTPUTS, INPUTS, f'python2 {SCRIPT} --run')
    NCR.MWF(ns, "chipseq_targets", OUTPUTS);

if 1:
    _ = '''
    #########################
    #### Make pileups  ######
    '''
    import src.util as _util
    _get_file = _util._get_file
    SCRIPT = "./src/plot_bigwig_pileup.py"
    FORCE = 1
    DICTS = [
        dict({
        "MODULE_FILE": "./src/plot_bigwig_pileup.py",
        "DATA":{
            "PEAK_FILE": ('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),

    #         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
            "OUTDIR": _util._get_output_file("OUTPUT/fig-2d"),
            'FORCE':FORCE,
            "BIGWIG_FILES":[
                _get_file('Mapped_data/192C/S17/ELF3myc-17C-ZT10_S17_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S18/ELF3myc-27C-ZT10_S18_Ath-TAIR10_RPKM.bw'),
                
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S17/ELF3myc-17C-ZT10_S17_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S18/ELF3myc-27C-ZT10_S18_Ath-TAIR10_RPKM.bw'
            ],
            "AXIS_DICTS":[
                {'set_ylim':[0,12],
                'set_ylabel':"log2( RPKM_at_peak )"},
                {'set_xlabel':"distance to peak(bp)",
                 'set_ylabel':"RPKM",
                 'set_ylim':[0,800],
                }
            ],
    }
    }),

    dict({
        "MODULE_FILE":"./src/plot_bigwig_pileup.py",
        "DATA":{
            "PEAK_FILE": ('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),
            
    #         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
            "OUTDIR": _util._get_output_file("OUTPUT/fig-2e"),
            'FORCE':FORCE,
            "BIGWIG_FILES":[
                _get_file('Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw'),
                
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw',
            ],

            "AXIS_DICTS":[
                {'set_ylim':[0,12],
                'set_ylabel':"log2( RPKM_at_peak )"},
                {'set_xlabel':"distance to peak(bp)",
                 'set_ylabel':"RPKM",
                 'set_ylim':[0,800],
                }
            ],
    }
    }),


    ### note ylim = [75,200]
    dict({
        "MODULE_FILE":"./src/plot_bigwig_pileup.py",
        "DATA":{
            "PEAK_FILE": ('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),
            
    #         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
            "OUTDIR": _util._get_output_file("OUTPUT/fig-2e-inlet"),
            'FORCE':FORCE,
            "BIGWIG_FILES":[
                _get_file('Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw'),
                _get_file('Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw'),
                
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S1/1487-17-ZT10_S1_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S2/1487-27-ZT10_S2_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S3/1488-17-ZT10_S3_Ath-TAIR10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/192C/S4/1488-27-ZT10_S4_Ath-TAIR10_RPKM.bw',
            ],

            "AXIS_DICTS":[
                {'set_ylim':[0,12],
                'set_ylabel':"log2( RPKM_at_peak )"},
                {'set_xlabel':"distance to peak(bp)",
                 'set_ylabel':"RPKM",
                 'set_ylim':[75,200],
                }
            ],
    }
    }),
        dict({
        "MODULE_FILE":"./src/plot_bigwig_pileup.py",
        "DATA":{
            "PEAK_FILE": ('OUTPUT/chipseq_differential_binding.peak_list.bed.summit'),
    #         "PEAK_FILE":'/home/feng/static/lists/1112__ELF3__chipTarg.narrowPeak',
            "OUTDIR": _util._get_output_file("OUTPUT/fig-2f"),
            'FORCE':FORCE,
            "BIGWIG_FILES":[
                _get_file('Mapped_data/189C/S10/gELF3myc-17C_S10_RPKM.bw'),
                _get_file('Mapped_data/189C/S11/gELF3myc-27C_S11_RPKM.bw'),
                _get_file('Mapped_data/189C/S16/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-17C_S16_RPKM.bw'),
                _get_file('Mapped_data/189C/S17/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-27C_S17_RPKM.bw'),                        
                
    #             _getFile("189C/S10/gELF3myc-17C_S10_RPKM.bw"),
               
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S10/gELF3myc-17C_S10_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S11/gELF3myc-27C_S11_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S16/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-17C_S16_RPKM.bw',
    #             '/home/feng/writable/teamkj/__backup/mapped-data/ChIP-seq/Mapped_data/189C/S17/1469-gELF3-myc-elf3-1xELF4-ox-ZT10-27C_S17_RPKM.bw'
            ],
                
            "AXIS_DICTS":[
                {'set_ylim':[0,12],
                'set_ylabel':"log2( RPKM_at_peak )"},
                {'set_xlabel':"distance to peak(bp)",
                 'set_ylabel':"RPKM",
                 'set_ylim':[0,800],
                }
            ],
    }
    }),
        
    ]    
    OUTS = []
    for DICT in DICTS[:]:
        OFNAME = f'{DICT["DATA"]["OUTDIR"]}/MAIN.png'
        #### DAG Spec
        RULE.M(ns, SCRIPT, None)
        [RULE.MWF(ns, x ) for x in DICT["DATA"]["BIGWIG_FILES"]]
        RULE.MWF(ns, 
            OFNAME,
            " ".join([SCRIPT, DICT["DATA"]["PEAK_FILE"]] + DICT["DATA"]["BIGWIG_FILES"]), 
            lambda c,DICT=DICT:LSC(f'''python2 -u - <<EOF\nfrom pymisca.atto_job import ModuleJob;ModuleJob({repr(DICT)})\nEOF'''))
        OUTS.append(OFNAME)

    NCR.M(ns, "chipseq_pileups", " ".join(OUTS))




if 1:
    ##### gene_lists
    CMD = f'''python2 -u - <<EOF
import pymisca.ext as pyext
from src.rnaseq_figure_2 import rnaseq_figure
import src.util as _util
df = pyext.OrderedDict()

template = rnaseq_figure
template.input_rnaseq="{input_rnaseq}" 
template.input_datasets_meta = "{input_datasets_meta}" 
template.input_markers_df = "{input_markers_df}" 
template.FNAME  = "NA0"
template.OFNAME = "NA1"
template.TARGET = "NA2"
job = rnaseq_figure()


df["signature_targets"] = (job.signature_targets)
df["chipseq_targets_genes"] = pyext.readData( "{{c.i[0]}}")['feat_acc'].unique()
df = pyext.pd.DataFrame.from_dict(df,orient='index').T
df.to_csv("{{c.o[0]}}", index =0 )
EOF
'''
    INPUTS = f"\
    OUTPUT/chipseq_targets_genes_job.peak_list.csv \
    src/rnaseq_figure_2.py src/util.py {input_rnaseq} {input_datasets_meta} {input_markers_df}"
    
    "pyext"
    OUTPUTS = _util._get_output_file("OUTPUT/gene_lists_dataframe.csv")
    RULE.MWF(ns,
        OUTPUTS,
        INPUTS,
        CMD
    )
    NCR.MWF(ns, "gene_lists", OUTPUTS)

NCR.M(ns, "OUTPUT", "boxplots heatmaps chipseq_targets chipseq_pileups gene_lists", )
NCR.MWF(ns, "clean", None, "rm -rf OUTPUT docs")


INPUT = "OUTPUT/"
OUTPUT = "docs/README.md"
CMD = f'''
DIR={INPUT}
TARGET=docs
mkdir -p $TARGET
mkdir -p _build_temp && mv -f chip* home* genes* infiles* job* temp* _temp* _build_temp || true


ln -f README.md $TARGET
cp -prf src/ -t $TARGET;
#rm -rf $TARGET/src/static/
cp -pfr {{c.i[1]}} -t $TARGET
cp -plrf $DIR -t $TARGET

find $TARGET -type l -delete  ### jekyll does not work with symlink
cd $TARGET
URL=https://gist.githubusercontent.com/glowinthedark/b1f5900be2490c5371f827a49fd09f49/raw/db45a596db2274b94206f4dfa6d479cee4e49845/generate_directory_index.py
curl -L $URL | python2 -

echo "[DONE]"
'''
RULE.MWF(ns, OUTPUT, "OUTPUT LUCKFILE.py.dot.svg", CMD)
NCR.MWF(ns, "docs", OUTPUT)

# NCR.MWF(ns, "_sink", "OUTPUT clean docs")



'REQUIRE_APT_GRAPHVIZ_DOT_BINARY'
from luck.graph import rules_to_graph
from graphviz import Digraph
# from graphviz import nohtml
import graphviz
import jinja2
from jinja2 import Template, StrictUndefined
import os
def jinja2_format(s,**context):
    # d = context.copy()
    d = __builtins__.copy()
    d.update(context)
    # .update(__builtins__)
    return Template(s,undefined=StrictUndefined).render(**d)

def rule_to_label( filename, ruletype, filesize):
    import os
    fmt = '''<       
        <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0">
        {#
          <TR>
            <TD ALIGN="LEFT" BGCOLOR="lightblue">attr</TD>
            <TD ALIGN="LEFT" BGCOLOR="lightblue">type</TD>
            <TD ALIGN="LEFT" BGCOLOR="lightblue">value</TD>
          </TR>
        #} 
          <TR>
            <TD ALIGN="LEFT" BGCOLOR="lightblue"> filename </TD>
            <TD ALIGN="LEFT" BGCOLOR="white" HREF="{{filename}}">
              <FONT COLOR="{% if os.path.exists(filename) %}blue{%else%}black{%endif%}"><U>{{filename}}</U></FONT></TD> 
          </TR>




          <TR>
            <TD ALIGN="LEFT" BGCOLOR="lightblue"> filesize </TD>
            <TD ALIGN="LEFT" BGCOLOR="white">{{filesize}}</TD>
          </TR>

          <TR>
            <TD ALIGN="LEFT" BGCOLOR="lightblue"> ruletype </TD>
            <TD ALIGN="LEFT" BGCOLOR="white">{{ruletype}}</TD>
          </TR>

        </TABLE>
    >'''
    return jinja2_format(fmt , **locals())



def rules_to_graph(rules, g, output_file, output_format):
    _ = '''
    Use RELATIVE path when plotting graph
    '''
    if g is None:
        g = Digraph('G', strict=True,
            # graph_attr=dict(
            # autosize="false", 
            # size="25.7,8.3!", resolution="100"
            # )
        )
        g.attr(rankdir='RL')    
    for v in rules:
        PWD = os.getcwd()
        for _output in v.output.split():
            _output = os.path.relpath( _output, PWD)
            filesize = os_stat_safe( _output).st_size
            if v.input:
                g.node( _output,  
                    label = rule_to_label(_output, type(v).__name__, filesize=filesize),
                    # level="middle_or_sink"),
                    shape='plaintext')      
            else:
                g.node( _output,  
                    label = rule_to_label(_output, type(v).__name__, filesize=filesize),
                    shape='plaintext')      
                continue
            for _input in v.input.split():
                _input  = os.path.relpath(_input,PWD)
                g.edge( _input, _output)
                # v.output, _input)
                # pprint([v.output,  rule_to_label(v)])
    if output_file is not None:
        res = g.render(filename = output_file, format = output_format)
        return res
    else:
        return g

from luck.header import os_stat_safe
OFNAME = "LUCKFILE.py.dot"
print(os.path.realpath(OFNAME+'.svg'))

RULE.M(
    ns,
    f'LUCKFILE.py.dot.svg',
    None,
    lambda c:rules_to_graph( ns.values(), None, OFNAME, "svg")
    )


# g.render()
    # v.input.split()])
# for v in ns.values():
    # pprint([v.output,  v.input.split()])
