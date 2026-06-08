# Configuration file for the Sphinx documentation builder.

from __future__ import annotations
import re
from pathlib import Path
# `brisc` is imported from site-packages of the active conda env (built
# from /home/wainberg/brisc/brisc/).  No sys.path mutation needed.

project = "brisc"
author = "Keon Arbabi & Michael Wainberg"
copyright = "2025, Keon Arbabi & Michael Wainberg"

templates_path = ["_templates"]

# -- General configuration ---------------------------------------------------

extensions = [
    "sphinx.ext.autodoc",
    "sphinx.ext.autosummary",
    "sphinx.ext.napoleon",
    "sphinx.ext.viewcode",
    "myst_parser",
    "sphinx_copybutton",
    "sphinx_design",
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
    "js/api-scrollspy.js",
]
html_show_sourcelink = False

html_theme_options = {
    "navbar_end": ["theme-switcher"],
    "show_version_warning_banner": False,
    # Drop the "On this page" secondary sidebar everywhere; a persistent
    # primary-sidebar toggle (added via JS/CSS) takes its place.
    "secondary_sidebar_items": [],
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

# -- Fix up Markdown-isms in docstrings on the fly -------------------------

# Match Markdown links in docstrings -- tolerates whitespace (including
# a line break) between `]` and `(`, and URLs without an http(s):// prefix.
_md_link_re = re.compile(
    r"\[(`?)([^\]]+?)\1\]\s*\(([^\s\)]+)\)",
    re.DOTALL,
)

# Markdown fenced code blocks (```lang ... ```) are not valid RST -- autodoc
# parses docstrings as reStructuredText, so the fence renders as literal
# "`lang" text. Convert them to RST literal blocks. (brisc's docstrings use
# Markdown fences throughout; fixing them here avoids editing the source.)
_fence_open_re = re.compile(r"^(\s*)```(\w*)\s*$")
_fence_close_re = re.compile(r"^\s*```\s*$")

def _md_fences_to_rst(lines):
    out = []
    i, n = 0, len(lines)
    while i < n:
        m = _fence_open_re.match(lines[i])
        if not m:
            out.append(lines[i])
            i += 1
            continue
        indent, lang = m.group(1), m.group(2)
        # Gather the fenced body up to the closing ``` (or end of docstring).
        body, j = [], i + 1
        while j < n and not _fence_close_re.match(lines[j]):
            body.append(lines[j])
            j += 1
        # Emit an RST code block, re-indenting the body 4 spaces past the
        # directive while preserving each line's relative indentation.
        out.append(f"{indent}.. code-block::{(' ' + lang) if lang else ''}")
        out.append("")
        for b in body:
            if not b.strip():
                out.append("")
            else:
                stripped = b[len(indent):] if b.startswith(indent) else b.lstrip()
                out.append(f"{indent}    {stripped}")
        out.append("")
        i = j + 1  # skip the closing fence
    return out

def _md_to_rst_links(app, what, name, obj, options, lines):
    # 1) Markdown code fences -> RST literal blocks.
    fenced = _md_fences_to_rst(lines)
    if fenced != lines:
        lines[:] = fenced
    # 2) Markdown links -> RST links.
    text = "\n".join(lines)
    def _repl(m):
        label, url = m.group(2), m.group(3)
        if not url.startswith(("http://", "https://")):
            url = "https://" + url
        return f"`{label} <{url}>`_"
    new_text = _md_link_re.sub(_repl, text)
    if new_text != text:
        lines[:] = new_text.split("\n")

# Build method → API doc URL mappings
def _build_api_links():
    """Build a dict of method_name → relative URL for all documented methods."""
    from brisc import SingleCell, Pseudobulk, DE
    links = {}
    # class names → their index pages
    links['SingleCell'] = 'api/single_cell/index.html'
    links['Pseudobulk'] = 'api/pseudobulk/index.html'
    links['DE'] = 'api/de/index.html'
    # methods and properties
    for cls, prefix in [(SingleCell, 'api/single_cell/api/brisc.SingleCell'),
                        (Pseudobulk, 'api/pseudobulk/api/brisc.Pseudobulk'),
                        (DE, 'api/de/api/brisc.DE')]:
        for name in dir(cls):
            if name.startswith('_') and name != '__init__':
                continue
            url = f'{prefix}.{name}.html'
            links[(cls.__name__, name)] = url
            # bare name → SingleCell takes priority, but don't overwrite class names
            if name not in links or cls is SingleCell:
                links[name] = url
    return links

# Built lazily on first use (see _get_api_links) to keep Sphinx startup fast —
# importing single_cell pulls in h5py/numpy/polars/pyarrow/scipy and takes ~2min.
_api_links = None

def _get_api_links():
    global _api_links
    if _api_links is None:
        _api_links = _build_api_links()
    return _api_links

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
    url = _get_api_links().get(name)
    if url:
        rel_url = '../' * depth + url
        return (f'<a href="{rel_url}" class="api-link" '
                f'style="text-decoration:none;color:inherit">'
                f'<span class="{css_class}">{name}</span></a>')
    return f'<span class="{css_class}">{name}</span>'

_SIDEBAR_TITLE_MAP = [
    ("api/single_cell/", "SingleCell"),
    ("api/pseudobulk/",  "Pseudobulk"),
    ("api/de/",          "DE"),
]
_SIDEBAR_TITLE_RE = re.compile(
    r'<p class="bd-links__title"[^>]*>Section Navigation</p>'
)
_SIDEBAR_ENTRY_RE = re.compile(
    r'>brisc\.(SingleCell|Pseudobulk|DE)\.([A-Za-z_][A-Za-z0-9_]*)<'
)

# Category basename (matches the rst filename) → Sphinx-generated heading
# anchor on the class's index page.
_CATEGORY_ANCHORS = {
    "constructor":           "constructor",
    "io":                    "i-o",
    "properties":            "properties",
    "data_access":           "data-access",
    "manipulation":          "manipulation",
    "structural":            "structural",
    "analysis":              "analysis",
    "utility":               "utility",
    "dictionary_interface":  "dictionary-interface",
    "transformation":        "transformation",
}
_CAT_NAMES = "|".join(re.escape(k) for k in _CATEGORY_ANCHORS)
# Matches sidebar toctree-l1 category entries so their href can be
# rewritten to point at the anchor on the class's index page.
_SIDEBAR_CAT_RE = re.compile(
    r'(<li class="toctree-l1[^"]*"><a class="reference internal" href=")'
    r'((?:\.\./)*)'
    rf'({_CAT_NAMES})\.html"'
)

def _semantic_highlight(app, exception=None):
    """Post-process HTML to add semantic classes and API links."""
    if exception is not None:
        return
    outdir = Path(app.builder.outdir)
    for html_file in outdir.rglob("*.html"):
        text = html_file.read_text()
        original = text
        rel = html_file.relative_to(outdir)
        depth = len(rel.parts) - 1  # e.g. tutorials/foo.html → depth 1
        rel_str = str(rel).replace("\\", "/")

        # -- Sidebar customization (runs on every page) --
        # Rename "Section Navigation" to the active class/module name.
        new_title = None
        for prefix, label in _SIDEBAR_TITLE_MAP:
            if rel_str.startswith(prefix):
                new_title = label
                break
        if new_title:
            text = _SIDEBAR_TITLE_RE.sub(
                f'<p class="bd-links__title" role="heading" aria-level="1">{new_title}</p>',
                text,
            )
        # Strip "single_cell.<Class>." prefix from sidebar entry labels
        # (and anywhere else that shows the fully-qualified dotted name).
        text = _SIDEBAR_ENTRY_RE.sub(r'>\2<', text)

        # Rewrite sidebar category links to jump to the anchored section on
        # the class index page (instead of loading a dedicated category
        # page).  e.g. href="constructor.html"  →  href="index.html#constructor"
        def _cat_repl(m):
            prefix, upward, cat = m.group(1), m.group(2), m.group(3)
            anchor = _CATEGORY_ANCHORS.get(cat, cat)
            return f'{prefix}{upward}index.html#{anchor}"'
        text = _SIDEBAR_CAT_RE.sub(_cat_repl, text)

        # API pages: mark the body so CSS can hide the right "On this page"
        # secondary sidebar -- the left-hand sidebar is the sole navigation
        # there, with category entries jumping to the anchored section on
        # the class index page (rewritten above).
        if rel_str.startswith("api/"):
            text = text.replace('<body ', '<body class="api-page" ', 1)

        # -- Code-block semantic highlighting (only where pygments ran) --
        if '<span class="n">' not in text:
            if text != original:
                html_file.write_text(text)
            continue

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

        # Scanpy-style params: split "name (type) – desc" into two lines.
        # Handles both the form with a description (em-dash present) and
        # the form without (just name + type).
        def _param_repl(m):
            name, types, sep = m.group(1), m.group(2), m.group(3)
            header = (
                f'<li><p class="param-header"><strong>{name}</strong> : '
                f'{types}</p>'
            )
            if '–' in sep:
                return header + '<p class="param-desc">'
            return header
        text = re.sub(
            r'<li><p><strong>([^<]+)</strong>\s*'
            r'\(([^)]*(?:\([^)]*\))*[^)]*)\)'
            r'(\s*–\s*|</p>)',
            _param_repl,
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

        # Collapse the Scalar-style union (str | int | float | Decimal | date
        # | time | datetime | timedelta | bool | bytes | Expr | Series | ...)
        # into the Scalar typedef link. Also handles the Iterable[...] variant.
        _scalar_inner = (
            r'<em>str</em><em> \| </em>'
            r'<em>int</em><em> \| </em>'
            r'<em>float</em><em> \| </em>'
            r'<em>Decimal</em><em> \| </em>'
            r'<em>date</em><em> \| </em>'
            r'<em>time</em><em> \| </em>'
            r'<em>datetime</em><em> \| </em>'
            r'<em>timedelta</em><em> \| </em>'
            r'<em>bool</em><em> \| </em>'
            r'<em>bytes</em><em> \| </em>'
            r'<em>Expr</em><em> \| </em>'
            r'<em>Series</em>'
        )
        _scalar_full = (
            _scalar_inner +
            r'(?:<em> \| </em><em>Iterable</em><em>\[</em>' +
            _scalar_inner +
            r'<em>\]</em>)?'
        )
        _scalar_link = (
            f'<a href="{"../" * depth}api/single_cell/typedefs.html"'
            f' style="text-decoration:none">'
            f'<em>Scalar</em></a> '
        )
        text = re.sub(_scalar_full, _scalar_link, text)

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
    app.connect("build-finished", _semantic_highlight, priority=901)
