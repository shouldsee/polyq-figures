
def SVG_CONTENT():
	from xml.dom import minidom
	svg_file = "LUCKFILE.py.dot.svg"
	doc = minidom.parse(svg_file)  # parseString also exists
	# path_strings = [path.getAttribute('d') for path
	#                 in doc.getElementsByTagName('path')]
	x = doc.getElementsByTagName('svg')[0]
	buf = ''.join([xx.toxml() for xx in x.childNodes])           
	# print(buf)
	return buf


import jinja2
from jinja2 import Template, StrictUndefined
def jinja2_format(s,**context):
    # d = context.copy()
    # d = __builtins__.copy()
    d = {}
    d.update(context)
    # .update(__builtins__)
    return Template(s,undefined=StrictUndefined).render(**d)

with open("index.html.jinja", 'r') as f:
	buf = jinja2_format(f.read(),**locals())
	print(buf)
