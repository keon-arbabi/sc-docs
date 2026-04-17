# Configuration file for the Sphinx documentation builder.

from __future__ import annotations
import re
import sys
from pathlib import Path
from markupsafe import Markup
sys.path.insert(0, '/home/karbabi')

project = "brisc"
author = "Keon Arbabi & Michael Wainberg"
copyright = Markup("2025,<br>Keon Arbabi &amp; Michael Wainberg")

templates_path = ["_templates"]

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
]

# Allow both .rst and .md
source_suffix = {
    ".rst": "restructuredtext",
    ".md": "markdown",
}

# MyST extensions for richer Markdown
myst_enable_extensions = [
    "colon_fence",
    "fieldlist",
    "deflist",
]

# Auto-generate anchors for H1-H3, so method-group headings are linkable
myst_heading_anchors = 3

# -- Autodoc / Autosummary configuration -------------------------------------

# Pull members in source order (matches the order in single_cell.py)
autodoc_member_order = "bysource"

# Move type hints from signature to parameter descriptions (Scanpy style)
autodoc_typehints = "description"

# Don't force one-parameter-per-line wrapping; let CSS handle natural wrapping
maximum_signature_line_length = 10000

# Type alias display is handled by the _semantic_highlight post-processor
autodoc_type_aliases = {}

# Generate stub files from autosummary directives
autosummary_generate = True

# Napoleon settings (Google-style docstrings with Args:, Returns:, Note:)
napoleon_google_docstring = True
napoleon_numpy_docstring = False
napoleon_include_init_with_doc = True
napoleon_include_private_with_doc = False
napoleon_use_admonition_for_notes = True
napoleon_use_param = True
napoleon_use_rtype = True
napoleon_preprocess_types = True

# -- Options for HTML output -------------------------------------------------

html_theme = "pydata_sphinx_theme"

html_static_path = ["_static"]
html_css_files = ["css/custom.css"]
html_js_files = [
    "js/benchmark-data.js",
    "js/benchmark.js",
    "js/carousel.js",
    "js/theme-fix.js",
]
html_show_sourcelink = False

html_theme_options = {
    "navbar_end": ["theme-switcher"],
    "show_version_warning_banner": False,
    "logo": {
        "image_dark": "_static/images/runner_logo_dark.svg",
        "image_light": "_static/images/runner_logo_light.svg",
        "text": "brisc documentation",
        "alt_text": "brisc documentation",
    },
}

# sidebar-nav-bs provides the full collapsible toctree navigation.
html_sidebars = {
    "**": ["sidebar-nav-bs"],
}

# Strip prompt prefixes from code copy
copybutton_prompt_text = r">>> |\.\.\. "
copybutton_prompt_is_regexp = True

# -- Convert Markdown links in docstrings to rST on the fly ----------------

_md_link_re = re.compile(r"\[(`?)([^\]]+?)\1\]\((https?://[^\)]+)\)")

def _md_to_rst_links(app, what, name, obj, options, lines):
    for i, line in enumerate(lines):
        lines[i] = _md_link_re.sub(r"`\2 <\3>`_", line)

def _overwrite_pygments_css(app, exception=None):
    """Replace pygments.css with Monokai Pro Ristretto Filter colors."""
    if exception is not None:
        return
    from _pygments_ristretto import RistrettoLightStyle, RistrettoDarkStyle
    from pygments.formatters import HtmlFormatter
    light = HtmlFormatter(style=RistrettoLightStyle)
    dark = HtmlFormatter(style=RistrettoDarkStyle)
    lines = []
    for fmt, prefix in [(light, 'html[data-theme="light"] .highlight'),
                        (dark,  'html[data-theme="dark"] .highlight')]:
        for line in fmt.get_linenos_style_defs():
            lines.append(f"{prefix} {line}")
        lines.extend(fmt.get_background_style_defs(prefix))
        lines.extend(fmt.get_token_style_defs(prefix))
    css_path = Path(app.builder.outdir) / "_static" / "pygments.css"
    css_path.write_text("\n".join(lines))

# Build method → API doc URL mappings
def _build_api_links():
    """Build a dict of method_name → relative URL for all documented methods."""
    from single_cell import SingleCell, Pseudobulk, DE
    links = {}
    # class names → their index pages
    links['SingleCell'] = 'api/single_cell/index.html'
    links['Pseudobulk'] = 'api/pseudobulk/index.html'
    links['DE'] = 'api/de/index.html'
    # methods and properties
    for cls, prefix in [(SingleCell, 'api/single_cell/api/single_cell.SingleCell'),
                        (Pseudobulk, 'api/pseudobulk/api/single_cell.Pseudobulk'),
                        (DE, 'api/de/api/single_cell.DE')]:
        for name in dir(cls):
            if name.startswith('_') and name != '__init__':
                continue
            url = f'{prefix}.{name}.html'
            links[(cls.__name__, name)] = url
            # bare name → SingleCell takes priority, but don't overwrite class names
            if name not in links or cls is SingleCell:
                links[name] = url
    return links

_api_links = _build_api_links()

# Pattern: <span class="n">NAME</span><span class="p">(...</span>  → function call
# The paren span may be `(`, `()`, or `(...)` depending on Pygments.
_call_re = re.compile(
    r'<span class="n">([^<]+)</span>'
    r'(<span class="p">\([^<]*</span>)')
# Pattern: <span class="n">NAME</span><span class="o">=</span>  → keyword arg
_kwarg_re = re.compile(
    r'<span class="n">([^<]+)</span>'
    r'(<span class="o">=</span>)')
# Pattern: .<span class="n">NAME</span><span class="p">(...</span>  → method call
_method_re = re.compile(
    r'(<span class="o">\.</span>)'
    r'<span class="n">([^<]+)</span>'
    r'(<span class="p">\([^<]*</span>)')

def _make_linked_span(name, css_class, depth):
    """Wrap a span in an <a> tag if the name is a known API method."""
    url = _api_links.get(name)
    if url:
        rel_url = '../' * depth + url
        return (f'<a href="{rel_url}" class="api-link" '
                f'style="text-decoration:none;color:inherit">'
                f'<span class="{css_class}">{name}</span></a>')
    return f'<span class="{css_class}">{name}</span>'

def _semantic_highlight(app, exception=None):
    """Post-process HTML to add semantic classes and API links."""
    if exception is not None:
        return
    outdir = Path(app.builder.outdir)
    for html_file in outdir.rglob("*.html"):
        text = html_file.read_text()
        if '<span class="n">' not in text:
            continue
        original = text
        # calculate relative depth for links
        rel = html_file.relative_to(outdir)
        depth = len(rel.parts) - 1  # e.g. tutorials/foo.html → depth 1

        # method calls: .name( → green + link
        def _method_repl(m):
            name = m.group(2)
            span = _make_linked_span(name, 'nf', depth)
            return m.group(1) + span + m.group(3)
        text = _method_re.sub(_method_repl, text)

        # function calls: name( → green + link
        def _call_repl(m):
            name = m.group(1)
            span = _make_linked_span(name, 'nf', depth)
            return span + m.group(2)
        text = _call_re.sub(_call_repl, text)

        # keyword args: name= → orange
        text = _kwarg_re.sub(
            r'<span class="na">\1</span>\2', text)

        # Scanpy-style params: split "name (type) – desc" into two lines
        text = re.sub(
            r'<li><p><strong>([^<]+)</strong>\s*'       # <li><p><strong>name</strong>
            r'\(([^)]*(?:\([^)]*\))*[^)]*)\)'           # (type, may have nested parens)
            r'\s*–\s*',                                  # –
            r'<li><p class="param-header"><strong>\1</strong> : \2</p>'
            r'<p class="param-desc">',
            text)

        # Simplify types: remove np.integer/np.floating/np.bool_ duplicates
        # "int | integer" → "int", "float | floating" → "float"
        for numpy_t, python_t in [('integer', 'int'), ('floating', 'float'),
                                   ('bool_', 'bool')]:
            # remove "| integer" or "integer | " patterns
            text = re.sub(
                rf'<em>\s*\|\s*</em><em>\s*{numpy_t}\s*</em>', '', text)
            text = re.sub(
                rf'<em>\s*{numpy_t}\s*</em><em>\s*\|\s*</em>', '', text)

        # Collapse SingleCellColumn union to typedef name
        _scc_pattern = (
            r'<em>str</em><em> \| </em><em>Expr</em><em> \| </em>'
            r'<em>Series</em><em> \| </em><em>ndarray</em><em> \| </em>'
            r'<em>Callable</em><em>\[</em><em>\[</em>'
            r'.*?SingleCell.*?'
            r'<em>\]</em><em>,\s*</em><em>Series</em><em> \| </em>'
            r'<em>ndarray</em><em>\]\s*</em>'
        )
        _scc_link = (
            f'<a href="{"../" * depth}api/single_cell/typedefs.html"'
            f' style="text-decoration:none">'
            f'<em>SingleCellColumn</em></a> '
        )
        text = re.sub(_scc_pattern, _scc_link, text)

        if text != original:
            html_file.write_text(text)

# -- Generate benchmark-data.js from sc-benchmarking CSVs ------------------

_BENCHMARK_DIR = Path("/home/karbabi/sc-benchmarking/output")
# For brisc's basic workflow, keep only PaCMAP to match scanpy/seurat which
# run a single embedding step.
_BASIC_BRISC_EXCLUDE = {
    "Embedding (LocalMAP)",
    "Embedding (UMAP)",
    "Embedding (UMAP hogwild)",
}

def _sum_timer_csv(path, exclude=None):
    import csv
    total = 0.0
    with open(path, newline="") as f:
        reader = csv.DictReader(f)
        for row in reader:
            op = row["operation"]
            if exclude and op in exclude:
                continue
            try:
                total += float(row["duration"])
            except (ValueError, TypeError):
                continue
    return total

def _generate_benchmark_data(app):
    import json
    cpu_workflows = [
        ("Basic workflow", "basic"),
        ("Label transfer", "transfer"),
        ("Pseudobulk differential expression", "de"),
    ]
    cpu_libs = [
        ("brisc",  "{prefix}_brisc_Parse_-1_timer.csv", True),
        ("scanpy", "{prefix}_scanpy_Parse_timer.csv",   False),
        ("seurat", "{prefix}_seurat_Parse_timer.csv",   False),
    ]
    groups = {}
    for label, prefix in cpu_workflows:
        bars = {}
        for lib_name, fmt, is_brisc in cpu_libs:
            csv_path = _BENCHMARK_DIR / fmt.format(prefix=prefix)
            if not csv_path.exists():
                continue
            exclude = _BASIC_BRISC_EXCLUDE if (is_brisc and prefix == "basic") else None
            bars[lib_name] = round(_sum_timer_csv(csv_path, exclude=exclude), 2)
        groups[label] = {"hardware": "cpu", "bars": bars}

    # GPU variant of the basic workflow: brisc vs rapids-single-cell on the
    # same 10M-cell Parse PBMC dataset but on 96 CPUs + 4x H100 GPUs.
    gpu_files = [
        ("brisc",  "basic_brisc_Parse_-1_gpu_timer.csv"),
        ("rapids", "basic_rapids_Parse_gpu_timer.csv"),
    ]
    gpu_bars = {}
    for lib_name, fname in gpu_files:
        csv_path = _BENCHMARK_DIR / fname
        if csv_path.exists():
            gpu_bars[lib_name] = round(_sum_timer_csv(csv_path), 2)
    if gpu_bars:
        groups["Basic workflow \u00b7 CPU vs GPU"] = {
            "hardware": "gpu",
            "note": "96 CPUs, 4\u00d7 H100 GPU, 752 GB RAM",
            "bars": gpu_bars,
        }

    payload = {
        "subtitle": "192 CPUs, 755 GB RAM",
        "groups": groups,
    }
    out_path = Path(app.srcdir) / "_static" / "js" / "benchmark-data.js"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(
        "window.BENCHMARK_DATA = " + json.dumps(payload, indent=2, ensure_ascii=False) + ";\n"
    )

def setup(app):
    app.connect("autodoc-process-docstring", _md_to_rst_links)
    app.connect("builder-inited", _generate_benchmark_data)
    app.connect("build-finished", _overwrite_pygments_css, priority=900)
    app.connect("build-finished", _semantic_highlight, priority=901)
