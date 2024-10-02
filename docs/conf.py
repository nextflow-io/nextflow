# -*- coding: utf-8 -*-
#
# Nextflow documentation build configuration file, created by
# sphinx-quickstart on Sat May  5 16:57:28 2012.
#
# This file is execfile()d with the current directory set to its containing dir.
#
# Note that not all possible configuration values are present in this
# autogenerated file.
#
# All configuration values have a default; values that are commented out
# serve to show the default.

import sys, os
import sphinx_rtd_theme
import datetime
import subprocess
import shlex

# If extensions (or modules to document with autodoc) are in another directory,
# add these directories to sys.path here. If the directory is relative to the
# documentation root, use os.path.abspath to make it absolute, like shown here.
#sys.path.insert(0, os.path.abspath('.'))

# -- General configuration -----------------------------------------------------

# If your documentation needs a minimal Sphinx version, state it here.
#needs_sphinx = '1.0'

# Add any Sphinx extension module names here, as strings. They can be extensions
# coming with Sphinx (named 'sphinx.ext.*') or your custom ones.
extensions = [
  'sphinx.ext.mathjax',
  'sphinxcontrib.mermaid',
  'sphinxext.rediraffe',
  'sphinx_rtd_theme',
  'myst_parser'
]

myst_enable_extensions = ['colon_fence', 'deflist', 'dollarmath']

myst_heading_anchors = 3

rediraffe_redirects = {
    'getstarted.md': 'install.md',
    'basic.md': 'overview.md',
    'tracing.md': 'reports.md',
    'mail.md': 'notifications.md',
    'operator.md': 'reference/operator.md',
    'dsl2.md': 'dsl1.md'
}

# Add any paths that contain templates here, relative to this directory.
templates_path = ['_templates']

# The suffix of source filenames.
source_suffix = '.md'

# The encoding of source files.
#source_encoding = 'utf-8-sig'

# The master toctree document.
master_doc = 'index'

# General information about the project.
project = u'Nextflow'
copyright = str(datetime.date.today().year) + u', Seqera Labs, S.L'

# The version info for the project you're documenting, acts as replacement for
# |version| and |release|, also used in various other places throughout the
# built documents.

# AUTO-GENERATED BY GIT
# Check if we are checked out at a specific version or not
process = subprocess.Popen(shlex.split("git tag --points-at HEAD"), stdout=subprocess.PIPE)

# The full version, including alpha/beta/rc tags.
release = process.communicate()[0].decode("utf-8").strip()

# The short X.Y version.
version = release.replace("-edge", "")

# The name for this set of Sphinx documents.  If None, it defaults to
# "<project> v<release> documentation".
html_title = f"Nextflow {release} documentation"

# Get the current sha if not checked out at a specific version
if len(release) == 0:
    process = subprocess.Popen(shlex.split("git rev-parse --short HEAD"), stdout=subprocess.PIPE)
    current_sha = process.communicate()[0].decode("utf-8").strip()
    release = f"dev <samp>({current_sha})</samp>"
    html_title = "Nextflow documentation"


# A shorter title for the navigation bar.  Default is the same as html_title.
#html_short_title = None

# The language for content autogenerated by Sphinx. Refer to documentation
# for a list of supported languages.
#language = None

# There are two options for replacing |today|: either, you set today to some
# non-false value, then it is used:
#today = ''
# Else, today_fmt is used as the format for a strftime call.
#today_fmt = '%B %d, %Y'

# List of patterns, relative to source directory, that match files and
# directories to ignore when looking for source files.
exclude_patterns = ['_build', '**README.md']

# The reST default role (used for this markup: `text`) to use for all documents.
#default_role = None

# If true, '()' will be appended to :func: etc. cross-reference text.
#add_function_parentheses = True

# If true, the current module name will be prepended to all description
# unit titles (such as .. function::).
#add_module_names = True

# If true, sectionauthor and moduleauthor directives will be shown in the
# output. They are ignored by default.
#show_authors = False

# The name of the Pygments (syntax highlighting) style to use.
pygments_style = 'default'

# A list of ignored prefixes for module index sorting.
#modindex_common_prefix = []


# -- Options for HTML output ---------------------------------------------------

# The theme to use for HTML and HTML Help pages.  See the documentation for
# a list of builtin themes.

html_theme = "sphinx_rtd_theme"

# Theme options are theme-specific and customize the look and feel of a theme
# further.  For a list of options available for each theme, see the
# documentation.
html_theme_options = {
  'logo_only': True,
  'display_version': False
}

html_context = {
    "display_github": True,
    "github_user": "nextflow-io",
    "github_repo": "nextflow",
    "github_version": "master",
    "conf_py_path": "/docs/"
}

# Nextflow theme
html_css_files = ['theme.css']

# Add any paths that contain custom themes here, relative to this directory.
#html_theme_path = []

# The name of an image file (relative to this directory) to place at the top
# of the sidebar.
html_logo = '_static/nextflow-logo.png'

# The name of an image file (within the static path) to use as favicon of the
# docs.  This file should be a Windows icon file (.ico) being 16x16 or 32x32
# pixels large.
html_favicon = '_static/favicon.ico'

# Add any paths that contain custom static files (such as style sheets) here,
# relative to this directory. They are copied after the builtin static files,
# so a file named "default.css" will overwrite the builtin "default.css".
html_static_path = ['_static']

# If not '', a 'Last updated on:' timestamp is inserted at every page bottom,
# using the given strftime format.
#html_last_updated_fmt = '%b %d, %Y'

# If true, SmartyPants will be used to convert quotes and dashes to
# typographically correct entities.
html_use_smartypants = False

# Custom sidebar templates, maps document names to template names.
html_sidebars = {
    '**': ['globaltoc.html']
}

# Additional templates that should be rendered to pages, maps page names to
# template names.
#html_additional_pages = {}

# If false, no module index is generated.
#html_domain_indices = True

# If false, no index is generated.
#html_use_index = True

# If true, the index is split into individual pages for each letter.
#html_split_index = False

# If true, links to the reST sources are added to the pages.
html_show_sourcelink = False

# If true, "Created using Sphinx" is shown in the HTML footer. Default is True.
#html_show_sphinx = True

# If true, "(C) Copyright ..." is shown in the HTML footer. Default is True.
#html_show_copyright = True

# If true, an OpenSearch description file will be output, and all pages will
# contain a <link> tag referring to it.  The value of this option must be the
# base URL from which the finished HTML is served.
#html_use_opensearch = ''

# This is the file name suffix for HTML files (e.g. ".xhtml").
#html_file_suffix = None

# Output file base name for HTML help builder.
htmlhelp_basename = 'Nextflowdoc'


# -- Options for LaTeX output --------------------------------------------------

latex_elements = {
# The paper size ('letterpaper' or 'a4paper').
#'papersize': 'letterpaper',

# The font size ('10pt', '11pt' or '12pt').
#'pointsize': '10pt',

# Additional stuff for the LaTeX preamble.
#'preamble': '',
}

# Grouping the document tree into LaTeX files. List of tuples
# (source start file, target name, title, author, documentclass [howto/manual]).
latex_documents = [
  ('index', 'Nextflow.tex', u'Nextflow Documentation',
   u'Paolo Di Tommaso', 'manual'),
]

# The name of an image file (relative to this directory) to place at the top of
# the title page.
latex_logo = '_static/nextflow-logo.png'

# For "manual" documents, if this is true, then toplevel headings are parts,
# not chapters.
#latex_use_parts = False

# If true, show page references after internal links.
#latex_show_pagerefs = False

# If true, show URL addresses after external links.
#latex_show_urls = False

# Documents to append as an appendix to all manuals.
#latex_appendices = []

# If false, no module index is generated.
#latex_domain_indices = True


# -- Options for manual page output --------------------------------------------

# One entry per manual page. List of tuples
# (source start file, name, description, authors, manual section).
man_pages = [
    ('index', 'nextflow', u'Nextflow Documentation',
     [u'Paolo Di Tommaso'], 1)
]

# If true, show URL addresses after external links.
#man_show_urls = False


# -- Options for Texinfo output ------------------------------------------------

# Grouping the document tree into Texinfo files. List of tuples
# (source start file, target name, title, author,
#  dir menu entry, description, category)
texinfo_documents = [
  ('index', 'Nextflow', u'Nextflow Documentation',
   u'Paolo Di Tommaso', 'Nextflow', 'One line description of project.',
   'Miscellaneous'),
]

# Documents to append as an appendix to all manuals.
#texinfo_appendices = []

# If false, no module index is generated.
#texinfo_domain_indices = True

# How to display URL addresses: 'footnote', 'no', or 'inline'.
#texinfo_show_urls = 'footnote'


# -- Options for Epub output ---------------------------------------------------

# Bibliographic Dublin Core info.
epub_title = u'Nextflow'
epub_author = u'Paolo Di Tommaso'
epub_publisher = u'Paolo Di Tommaso'
epub_copyright = u'Copyright 2013-2024, Seqera Labs'

# The language of the text. It defaults to the language option
# or en if the language is not set.
#epub_language = ''

# The scheme of the identifier. Typical schemes are ISBN or URL.
#epub_scheme = ''

# The unique identifier of the text. This can be a ISBN number
# or the project homepage.
#epub_identifier = ''

# A unique identification for the text.
#epub_uid = ''

# A tuple containing the cover image and cover page html template filenames.
#epub_cover = ()

# HTML files that should be inserted before the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_pre_files = []

# HTML files that should be inserted after the pages created by sphinx.
# The format is a list of tuples containing the path and title.
#epub_post_files = []

# A list of files that should not be packed into the epub file.
#epub_exclude_files = []

# The depth of the table of contents in toc.ncx.
#epub_tocdepth = 3

# Allow duplicate toc entries.
#epub_tocdup = True

# see also: https://github.com/pygments/pygments/blob/master/pygments/lexers/jvm.py
from pygments.lexer import RegexLexer
from sphinx.highlighting import lexers

class NextflowLexer(RegexLexer):
    """
    For Nextflow source code.
    """

    import re
    from pygments.lexer import bygroups, using, this, default
    from pygments.token import Comment, Operator, Keyword, Name, String, Number, Whitespace
    from pygments.util import shebang_matches

    name = 'Nextflow'
    url = 'https://nextflow.io/'
    aliases = ['nextflow']
    filenames = ['*.nf']
    mimetypes = ['text/x-nextflow']
    # version_added = '1.5'

    flags = re.MULTILINE | re.DOTALL

    tokens = {
        'root': [
            # Nextflow allows a file to start with a shebang
            (r'#!(.*?)$', Comment.Preproc, 'base'),
            default('base'),
        ],
        'base': [
            (r'[^\S\n]+', Whitespace),
            (r'(//.*?)(\n)', bygroups(Comment.Single, Whitespace)),
            (r'/\*.*?\*/', Comment.Multiline),
            # keywords: go before method names to avoid lexing "throw new XYZ"
            # as a method signature
            (r'(assert|catch|else|'
             r'if|instanceof|new|return|throw|try|in|as)\b',
             Keyword),
            # method names
            (r'^(\s*(?:[a-zA-Z_][\w.\[\]]*\s+)+?)'  # return arguments
             r'('
             r'[a-zA-Z_]\w*'                        # method name
             r'|"(?:\\\\|\\[^\\]|[^"\\])*"'         # or double-quoted method name
             r"|'(?:\\\\|\\[^\\]|[^'\\])*'"         # or single-quoted method name
             r')'
             r'(\s*)(\()',                          # signature start
             bygroups(using(this), Name.Function, Whitespace, Operator)),
            (r'@[a-zA-Z_][\w.]*', Name.Decorator),
            (r'(def|enum|include|from|output|process|workflow)\b', Keyword.Declaration),
            (r'(boolean|byte|char|double|float|int|long|short|void)\b',
             Keyword.Type),
            (r'(true|false|null)\b', Keyword.Constant),
            (r'""".*?"""', String.Double),
            (r"'''.*?'''", String.Single),
            (r'"(\\\\|\\[^\\]|[^"\\])*"', String.Double),
            (r"'(\\\\|\\[^\\]|[^'\\])*'", String.Single),
            (r'/(\\\\|\\[^\\]|[^/\\])*/', String),
            (r"'\\.'|'[^\\]'|'\\u[0-9a-fA-F]{4}'", String.Char),
            (r'(\.)([a-zA-Z_]\w*)', bygroups(Operator, Name.Attribute)),
            (r'[a-zA-Z_]\w*:', Name.Label),
            (r'[a-zA-Z_$]\w*', Name),
            (r'[~^*!%&\[\](){}<>|+=:;,./?-]', Operator),
            (r'[0-9][0-9]*\.[0-9]+([eE][0-9]+)?[fd]?', Number.Float),
            (r'0x[0-9a-fA-F]+', Number.Hex),
            (r'[0-9]+L?', Number.Integer),
            (r'\n', Whitespace)
        ],
    }

    def analyse_text(text):
        return shebang_matches(text, r'nextflow')

lexers['nextflow'] = NextflowLexer()
