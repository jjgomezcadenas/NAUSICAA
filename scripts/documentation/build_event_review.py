#!/usr/bin/env python3
"""Build the offline event-analysis design, Petit appendices and short roadmap.

Python 3.11+, standard library only. Rebuild with make event-review.
Use --petit-root to select another local Petit checkout.
"""

from pathlib import Path as _Path
import sys as _sys
_REPO_ROOT = _Path(__file__).resolve().parents[2]
_sys.path[:0] = [str(_REPO_ROOT), str(_REPO_ROOT / "src")]

import argparse
from datetime import datetime, timezone
from hashlib import sha256
from html import escape
from pathlib import Path
import re
import subprocess

from docs.event_analysis_review.content.petit_review_content import assessment_pages
from docs.event_analysis_review.content.event_design_content import design_pages, roadmap_markdown

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / 'docs/event_analysis_review/index.html'


def code(value):
    return '<pre><code>' + escape(value.strip()) + '</code></pre>'


def table(headers, rows):
    return '<div class="table-scroll"><table><thead><tr>' + ''.join(
        f'<th scope="col">{h}</th>' for h in headers
    ) + '</tr></thead><tbody>' + ''.join(
        '<tr>' + ''.join(f'<td>{v}</td>' for v in row) + '</tr>' for row in rows
    ) + '</tbody></table></div>'


def section(title, content, anchor=None):
    attr = f' id="{anchor}"' if anchor else ''
    return f'<section class="card"{attr}><h2 tabindex="-1">{title}</h2>{content}</section>'


class Sources:
    def __init__(self, root):
        self.root = root
        self.files = {}

    def read(self, relative):
        if relative not in self.files:
            self.files[relative] = (self.root / relative).read_text()
        return self.files[relative]

    def excerpt(self, relative, marker, count=30):
        lines = self.read(relative).splitlines()
        first = next((i for i, line in enumerate(lines) if line.startswith(marker)), None)
        if first is None:
            raise ValueError(f'Missing source marker {marker!r} in {relative}')
        selected = lines[first:first + count]
        numbered = '\n'.join(f'{i + first + 1:4d}  {line}' for i, line in enumerate(selected))
        label = f'{relative}:{first + 1}–{first + len(selected)}'
        return '<details class="source"><summary>Source excerpt · ' + escape(label) + '</summary>' + code(numbered) + '</details>'

    def declaration(self, relative, name):
        source = self.read(relative)
        match = re.search(rf'^struct {name}\n.*?^end$', source, re.M | re.S)
        if not match:
            raise ValueError(f'Missing struct {name} in {relative}')
        return code(match.group()) + self.excerpt(relative, f'struct {name}', len(match.group().splitlines()))

    def provenance(self):
        result = subprocess.run(['git', '-C', str(self.root), 'rev-parse', 'HEAD'],
                                capture_output=True, text=True)
        commit = result.stdout.strip() if result.returncode == 0 else 'unavailable'
        rows = [(escape(path), f'<code>{sha256(text.encode()).hexdigest()}</code>')
                for path, text in sorted(self.files.items())]
        return ('<p>Source checkout: <code>' + escape(str(self.root)) + '</code>. '
                'Git HEAD: <code>' + escape(commit) + '</code>. '
                'Excerpts and file hashes describe the working-tree contents inspected for this page.</p>'
                '<details><summary>Inspected files and SHA-256 hashes</summary>'
                + table(['Petit source file', 'SHA-256'], rows) + '</details>')


def build():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--petit-root', type=Path, default=ROOT.parent / 'Petit')
    args = parser.parse_args()
    sources = Sources(args.petit_root.resolve())
    appendices = assessment_pages(sources, section, table, code)
    for page, title in zip(appendices, ('A. Petit implementation', 'B. Petit interfaces', 'C. Decision record')):
        page['title'] = title
    appendices[-1]['summary'] = 'Approved adaptation decisions and the validation proposal, with source context.'
    appendices[-1]['status'] = 'Decision record · analysis steps 1–4 implemented'
    pages = design_pages(ROOT, section, table, code) + appendices
    license_text = sources.read('LICENSE')
    provenance = sources.provenance()
    provenance += '<details><summary>Petit source licence · MIT</summary>' + code(license_text) + '</details>'
    local_files = ['config/event_analysis.toml', 'config/detector.toml', 'config/ion_tracks.toml', 'Project.toml',
                   'src/event_analysis/types.jl', 'src/event_analysis/coordinates.jl',
                   'src/ion_signal/event_io.jl', 'docs/event_analysis_review/content/event_design_content.py',
                   'docs/event_analysis_review/content/petit_review_content.py', 'scripts/documentation/build_event_review.py',
                   'docs/event_analysis_review/content/event_viewer_content.py', 'config/event_viewer.toml',
                   'docs/event_analysis_review/template.html']
    local_files += [str(p.relative_to(ROOT)) for p in sorted((ROOT / 'src/event_analysis').rglob('*')) if p.is_file() and str(p.relative_to(ROOT)) not in local_files]
    provenance += '<details><summary>NAUSICAA design inputs and documentation sources</summary>' + table(
        ['Local source', 'SHA-256'], [(escape(p), '<code>' + sha256((ROOT / p).read_bytes()).hexdigest() + '</code>') for p in local_files]) + '</details>'
    provenance += '<p class="caption">Built ' + datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC') + '. Thermal examples derive from the detector card. Steps one through four are linked in the checkpoints. Paired scans expose histogram sensitivity; production accuracy and the full-sample run remain for review.</p>'
    panels = []
    primary_links = []
    topic_groups = []
    for i, page in enumerate(pages):
        body = page['body']
        if page['id'] == 'active':
            body += section('Source provenance', provenance, 'provenance')
        sections = re.findall(r'<section class="card" id="([^"]+)"><h2[^>]*>(.*?)</h2>', body)
        toc = '<nav class="local-toc" aria-label="On this page"><strong>On this page</strong><ul>' + ''.join(
            f'<li><a href="#{anchor}">{title}</a></li>' for anchor, title in sections) + '</ul></nav>'
        primary_links.append(f'<li><a href="#{page["id"]}" data-search="{escape((page["title"] + " " + page["summary"]).lower(), quote=True)}">{page["title"]}</a></li>')
        topic_links = ''.join(f'<li><a href="#{anchor}" data-search="{escape((title + " " + page["title"]).lower(), quote=True)}">{title}</a></li>' for anchor, title in sections)
        topic_groups.append(f'<details class="topic-group"><summary>{page["title"]}</summary><ul class="review-links compact-links">{topic_links}</ul></details>')
        nav = '<nav class="page-nav" aria-label="Review sequence">'
        if i:
            nav += f'<a href="#{pages[i-1]["id"]}">← {pages[i-1]["title"]}</a>'
        if i + 1 < len(pages):
            nav += f'<a class="next" href="#{pages[i+1]["id"]}">Next: {pages[i+1]["title"]} →</a>'
        nav += '</nav>'
        panels.append(f'<article id="{page["id"]}" class="review-page" aria-labelledby="{page["id"]}-title">'
                      f'<header class="page-header"><div class="eyebrow">Design guide · {i+1} / {len(pages)}</div>'
                      f'<h1 id="{page["id"]}-title" tabindex="-1">{page["title"]}</h1>'
                      f'<p class="lead">{page["summary"]}</p><span class="badge">{page["status"]}</span></header>{toc}{body}{nav}</article>')
    template = (OUT.parent / 'template.html').read_text()
    assert template.count('<!-- REVIEW_PAGES -->') == 1
    for marker in ('<!-- PAGE_LINKS -->', '<!-- TOPIC_LINKS -->'):
        assert template.count(marker) == 1
    OUT.write_text(template.replace('<!-- REVIEW_PAGES -->', '\n'.join(panels))
                   .replace('<!-- PAGE_LINKS -->', '\n'.join(primary_links))
                   .replace('<!-- TOPIC_LINKS -->', '\n'.join(topic_groups)))
    (ROOT / 'docs/event_analysis_plan.md').write_text(roadmap_markdown())
    print(f'Built {OUT.relative_to(ROOT)}: {len(pages)} chapters; {len(sources.files)} inspected Petit files; roadmap refreshed.')


if __name__ == '__main__':
    build()
