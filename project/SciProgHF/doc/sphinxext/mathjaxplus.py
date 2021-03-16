# coding=utf-8

"""
This extensions adds ability to download math equations in latex format.
"""

import os

from docutils import nodes

try:
    from hashlib import sha1 as sha
except ImportError:
    from sha import sha

from sphinx.util.osutil import ensuredir
from sphinx.application import ExtensionError
from sphinx.ext.mathbase import setup_math as mathbase_setup

LATEX_TEMPLATE = r'''
\begin{equation}
%s
\end{equation}
'''


def get_latex_filename(content):
    """
    Returns filename for given content using sha1 hash of the content.
    """
    return '%s.tex' % sha(content).hexdigest()


def write_latex_file(file_path, content):
    """
    Writes content wrapped in latex equation to given file.
    """
    fp = open(file_path, 'w')
    fp.write(LATEX_TEMPLATE.strip() % content)
    fp.close()


def generate_latex_file(self, node):
    """
    Generates latex file for given math node and returns its path.
    """
    content = node['latex'].encode('utf-8')
    filename = get_latex_filename(content)
    file_path = os.path.join(self.builder.outdir, '_downloads', 'math', filename)
    file_rel_path = os.path.join(self.builder.dlpath, 'math', filename)

    ensuredir(os.path.dirname(file_path))
    write_latex_file(file_path, content)

    return file_rel_path


def html_visit_math(self, node):
    """
    Node visitor for math environment.
    """
    self.body.append(self.starttag(node, 'span', '', CLASS='math'))
    self.body.append(self.builder.config.mathjax_inline[0] +
                     self.encode(node['latex']) +
                     self.builder.config.mathjax_inline[1] + '</span>')
    raise nodes.SkipNode


def html_visit_displaymath(self, node):
    """
    Node visitor for displaymath environment.
    """
    self.body.append(self.starttag(node, 'div', CLASS='math'))
    if node['nowrap']:
        self.body.append(self.builder.config.mathjax_display[0] +
                         node['latex'] +
                         self.builder.config.mathjax_display[1])
        self.body.append('</div>')
        raise nodes.SkipNode

    parts = [prt for prt in node['latex'].split('\n\n') if prt.strip()]
    for i, part in enumerate(parts):
        part = self.encode(part)
        if i == 0:
            # necessary to e.g. set the id property correctly
            if node['number']:
                self.body.append('<span class="eqno">(%s)</span>' %
                                 node['number'])
        if '&' in part or '\\\\' in part:
            self.body.append(self.builder.config.mathjax_display[0] +
                             '\\begin{split}' + part + '\\end{split}' +
                             self.builder.config.mathjax_display[1])
        else:
            self.body.append(self.builder.config.mathjax_display[0] + part +
                             self.builder.config.mathjax_display[1])

    if self.builder.config.mathjax_generate_latex:
        latex = generate_latex_file(self, node)
        if latex:
            title = self.builder.config.mathjax_latex_title
            style = ''
            if self.builder.config.mathjax_latex_style:
                style = ' style="%s"' % self.builder.config.mathjax_latex_style
            self.body.append('<div%s><a href="%s">%s</a></div>\n' % (style, latex, title))
    self.body.append('</div>\n')

    raise nodes.SkipNode


def builder_inited(app):
    """
    Adds mathjax javascript file to app.
    """
    if not app.config.mathjax_path:
        raise ExtensionError('mathjax_path config value must be set for the '
                             'mathjax extension to work')
    app.add_javascript(app.config.mathjax_path)


def setup(app):
    mathbase_setup(app, (html_visit_math, None), (html_visit_displaymath, None))
    app.add_config_value('mathjax_path',
                         'http://cdn.mathjax.org/mathjax/latest/MathJax.js?'
                         'config=TeX-AMS-MML_HTMLorMML', False)
    app.add_config_value('mathjax_inline', [r'\(', r'\)'], 'html')
    app.add_config_value('mathjax_display', [r'\[', r'\]'], 'html')
    app.add_config_value('mathjax_generate_latex', True, 'html')
    app.add_config_value('mathjax_latex_title', 'Download as Latex', 'html')
    app.add_config_value('mathjax_latex_style', 'text-align:right;font-size:14px;', 'html')
    app.connect('builder-inited', builder_inited)
