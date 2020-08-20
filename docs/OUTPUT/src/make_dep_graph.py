from graphviz import Digraph
import pymisca.ext as pyext
from pymisca.ext import jf2 as _jf2
import json
from pymisca.atto_jobs import Shellexec
# from src.util import 
import src.util as _util

def graph_add_node_dict(s,name,d,**kw):
    if hasattr(d,'items',):
        d = d.items()
    return graph_add_node_table(s, name, d,**kw)


def graph_add_node_table(s, name, tab,**kw):
    name = str(name)
#     fmt = lambda x:str(x).replace('/','.')
    fmt = lambda x:x
    buf = _jf2(
           '''<
    <TABLE BORDER="0" CELLBORDER="1" CELLSPACING="0" CELLPADDING="4">
      {% for row in tab %}<TR>{% for value in row %}
      <TD PORT="{{fmt(value)}}" >{{(value)}}</TD>{% endfor %}  </TR>{% endfor %}
    </TABLE>>''')
    s.node( name,  buf,**kw)
    return s

# FNAME = 'OUTPUT/index.json.list'
def make_dep_graph(FNAME=_util.INDEX_FILE, save=1):
    DOT_FILE = FNAME+'.dot'
    FORMAT='svg'
    OFNAME= _util.get_output_file(DOT_FILE+'.'+FORMAT)
    it = [ json.loads(x.rstrip(),object_pairs_hook=pyext._DICT_CLASS) for x in open(FNAME,'r')]
    s = Digraph('structs', node_attr={'shape': 'plaintext'},
                graph_attr={ "rankdir":"LR" }
               )

    with s.subgraph(name='left') as s_left:
        with s.subgraph(name='right') as s_right:
            for i,d in enumerate(it):
            #     graph_add_node_dict(s,i,d)
                graph_add_node_dict(s_left,
                                    d['OUTPUT_FILE'],
                                    {'OUTPUT_FILE':d['OUTPUT_FILE']},
                                   href=d['OUTPUT_FILE'])
                graph_add_node_dict(s_right,
                                    d['RUNTIME_FILE'],
                                    {'RUNTIME_FILE':d['RUNTIME_FILE']},
                                   href=d['RUNTIME_FILE'])
                s.edge(d['OUTPUT_FILE'],d['RUNTIME_FILE'])
    s.save(DOT_FILE)
    res = Shellexec({'CMD_LIST':['dot',
                                 '-T%s'%FORMAT, DOT_FILE,'>',
                                 OFNAME,
                                ],'OUTDIR':'.'})       
     
    return s
if __name__ == '__main__':
    make_dep_graph()

#     FILE = s.render(format='svg')

# if __name__ == '__main__':
    if 0:
        s = Digraph('structs', node_attr={'shape': 'plaintext'})
        graph_add_node_table(s,'struct1', {"f1.1":"int", "f2":"value","f3":"value"}.items())
        graph_add_node_dict(s, 'struct2', {"f0":"one","f2":2})
        graph_add_node_dict(s, 'struct3', {"f0":"one","f2":2})

        s.edges([('struct1:f1', 'struct2:f0'), 
                 ('struct1:f2', 'struct3:f2')])

        FILE = s.render(format='svg',filename='test.svg')
        pyext.ipd.display(pyext.ipd.SVG(FILE))
    
    