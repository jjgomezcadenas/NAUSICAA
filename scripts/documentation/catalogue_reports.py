#!/usr/bin/env python3
"""Organize generated reports and build their searchable local HTML catalogue.

The registry records scientific subjects and verified workflow commands. Document
discovery requires a TeX documentclass, so figure PDFs and included tables are
not counted as reports. Existing provenance files are read but never rewritten.
Recorded report paths are resolved through the registry at read boundaries.
"""
from __future__ import annotations

import argparse
import hashlib
import html
import json
import os
from pathlib import Path
import re
import shutil

ROOT = Path(__file__).resolve().parents[2]
OUT = ROOT / "latex_generated"
REGISTRY = ROOT / "docs/generated_reports.json"
GROUPS = {
    "detector": "Detector and ELM",
    "transport": "Carrier transport",
    "reconstruction": "Event reconstruction",
    "ion_readout": "Ion readout and CSA",
    "electron_readout": "Electron readout and SiPM",
    "production": "Production and combined response",
    "literature": "Literature reviews",
    "archive": "Historical reference outputs",
}
INTERMEDIATES = {".aux", ".log", ".out", ".fls", ".fdb_latexmk", ".blg", ".bbl", ".toc", ".stdout"}


def relocate(registry):
    """Move output directories with byte-for-byte preservation."""
    moved = []
    for slug, spec in registry.items():
        old = OUT / slug
        new = OUT / spec["group"] / slug
        if old.is_symlink():
            if old.resolve() != new.resolve():
                raise ValueError(f"Unexpected compatibility link: {old}")
            old.unlink()
            continue
        if not old.exists():
            continue
        if new.exists():
            raise FileExistsError(new)
        before = {str(p.relative_to(old)): hashlib.sha256(p.read_bytes()).hexdigest()
                  for p in old.rglob("*") if p.is_file()}
        new.parent.mkdir(parents=True, exist_ok=True)
        old.rename(new)
        after = {str(p.relative_to(new)): hashlib.sha256(p.read_bytes()).hexdigest()
                 for p in new.rglob("*") if p.is_file()}
        if before != after:
            raise RuntimeError(f"Content changed while moving {slug}")
        moved.append(slug)
    return moved


def clean_intermediates(registry):
    """Keep compilation products under work/latex, preserving unique files."""
    count = 0
    for slug, spec in registry.items():
        directory = OUT / spec["group"] / slug
        if not directory.exists():
            continue
        for path in directory.rglob("*"):
            if not path.is_file() or not (path.suffix in INTERMEDIATES or
                                          path.name.endswith(".synctex.gz")):
                continue
            target = ROOT / "work/latex" / spec["group"] / slug / path.relative_to(directory)
            target.parent.mkdir(parents=True, exist_ok=True)
            if target.exists() and target.read_bytes() != path.read_bytes():
                digest = hashlib.sha256(path.read_bytes()).hexdigest()[:10]
                target = target.with_name(target.name + "." + digest)
            shutil.move(path, target)
            count += 1
    return count


def plain_tex(text):
    text = re.sub(r"\\(?:texttt|textbf|emph|mathrm|textrm|detokenize)\{([^{}]*)\}", r"\1", text)
    text = text.replace(r"\_", "_").replace(r"\%", "%").replace("~", " ")
    text = re.sub(r"\\[a-zA-Z]+\*?(?:\[[^]]*\])?", "", text)
    return re.sub(r"\s+", " ", text.replace("{", "").replace("}", "")).strip()


def brace_argument(text, command):
    match = re.search(re.escape(command) + r"\s*\{", text)
    if not match:
        return ""
    start = match.end()
    depth = 1
    for index in range(start, len(text)):
        if text[index] == "{": depth += 1
        elif text[index] == "}": depth -= 1
        if depth == 0: return plain_tex(text[start:index])
    return ""


def generate(registry):
    links = []
    def link(path, label):
        path = Path(path)
        if not path.is_absolute(): path = ROOT / path
        if not path.exists(): return ""
        links.append(path)
        return f'<a href="{html.escape(os.path.relpath(path, OUT))}">{html.escape(label)}</a>'

    scripts = list((ROOT / "scripts").rglob("*.py")) + list((ROOT / "scripts").rglob("*.jl"))
    script_text = {p: p.read_text(errors="replace") for p in scripts}
    entries = []
    for slug, spec in registry.items():
        directory = OUT / spec["group"] / slug
        if not directory.exists(): directory = OUT / slug
        if not directory.exists(): continue
        documents = []
        texts = []
        for tex in sorted(directory.rglob("*.tex")):
            if tex.is_symlink(): continue
            text = tex.read_text(errors="replace")
            if r"\documentclass" not in text: continue
            title = brace_argument(text, r"\title") or tex.stem.replace("_", " ")
            documents.append({"title": title, "tex": str(tex.relative_to(ROOT)),
                              "pdf": str(tex.with_suffix(".pdf").relative_to(ROOT)) if tex.with_suffix(".pdf").exists() else None})
            texts.append(text)
        for source_name in ([spec["source"]] if spec.get("source") else []) + spec.get("additional_sources", []):
            source = ROOT / source_name
            pdf = directory / source.with_suffix(".pdf").name
            if source.exists() and pdf.exists() and not any(d["pdf"] == str(pdf.relative_to(ROOT)) for d in documents):
                text = source.read_text(errors="replace")
                documents.append({"title": brace_argument(text, r"\title") or spec.get("title", source.stem),
                                  "tex": str(source.relative_to(ROOT)), "pdf": str(pdf.relative_to(ROOT))})
                texts.append(text)
        generators = [str(p.relative_to(ROOT)) for p,t in script_text.items()
                      if "latex_generated/" + slug in t or "latex_generated/" + spec["group"] + "/" + slug in t
                      or ("latex_generated" in t and re.search(r"['\"]"+re.escape(slug)+r"['\"]", t))]
        provenance = [str(p.relative_to(ROOT)) for p in sorted(directory.rglob("*"))
                      if p.is_file() and p.name in {"sources.json", "provenance.json", "sources.txt", "manifest.toml"}]
        # Link concrete local references found in generated documents. These are
        # evidence, not inferred physical settings or reconstructed run cards.
        refs = sorted(set(re.findall(r"(?:results|config|docs|notes)/[A-Za-z0-9_./-]+\.(?:h5|toml|json|tex)", "\n".join(texts))))
        refs = [p for p in refs if (ROOT/p).exists()]
        entries.append(dict(slug=slug, **spec, documents=documents, generators=generators,
                            provenance=provenance, references=refs))
    parts = ['<!doctype html><html lang="en"><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">',
             '<title>NAUSICAA · Generated notes</title><style>body{font:16px/1.55 system-ui;margin:0;color:#203040;background:#f6f7f8}main{max-width:1180px;margin:40px auto;padding:0 24px}h1{margin-bottom:0}h2{padding-top:20px;border-top:2px solid #dbe2e8}article{background:white;border:1px solid #dae1e6;border-radius:8px;padding:18px 22px;margin:14px 0}h3{margin:0}a{color:#126687}a+a{margin-left:14px}small,.muted{color:#586977}code{white-space:pre-wrap;background:#f0f3f5;padding:2px 5px}input{font:inherit;padding:10px;width:95%;margin:18px 0}nav a{display:inline-block;margin:3px 12px 3px 0}summary{cursor:pointer}li{margin:6px 0}.tag{font-size:12px;border-radius:4px;background:#e9f0f5;padding:3px 7px;float:right}.doc{margin:8px 0}footer{margin:35px 0}</style><main>',
             '<h1>NAUSICAA · Generated notes</h1><p>Reports, model notes and validation studies, grouped by scientific subject. Each entry links its saved documents and available provenance.</p>',
             '<nav>'+''.join(f'<a href="#{g}">{label}</a>' for g,label in GROUPS.items())+'</nav>',
             '<input id="search" type="search" placeholder="Filter by title, configuration, species or workflow…" aria-label="Filter reports">']
    for group, label in GROUPS.items():
        rows = [r for r in entries if r["group"] == group]
        if not rows: continue
        parts.append(f'<section id="{group}"><h2>{label} <small>({len(rows)})</small></h2>')
        for row in rows:
            title = row.get("title") or (row["documents"][0]["title"] if row["documents"] else row["slug"].replace("_", " ").capitalize())
            parts.append(f'<article><span class="tag">{html.escape(row["status"])}</span><h3>{html.escape(title)}</h3>')
            if row.get("description"): parts.append('<p>'+html.escape(row["description"])+'</p>')
            if row.get("configuration"): parts.append('<p><b>Configuration:</b> '+html.escape(row["configuration"])+'</p>')
            for document in row["documents"]:
                parts.append('<div class="doc">'+html.escape(document['title'])+' · '+link(document['tex'],'TeX')+' '+(link(document['pdf'],'PDF') if document['pdf'] else '<small>TeX available</small>')+'</div>')
            if row.get("html"): parts.append(link(row["html"],"Interactive HTML report"))
            extra = [link(row[k], k.capitalize()) for k in ("source", "result") if row.get(k)]
            extra += [link(p, Path(p).name) for p in row.get("additional_sources", [])]
            extra += [link(p,"Card: "+Path(p).name) for p in row.get("cards",[])]
            parts.append('<p>'+" ".join(extra)+'</p>')
            for key,label2 in (("presentation","Presentation only"),("workflow","Complete workflow")):
                if row.get(key): parts.append(f'<p>{label2}: <code>{html.escape(row[key])}</code></p>')
            if row.get("status_basis"): parts.append('<p>'+html.escape(row['status_basis'])+'</p>')
            details = row['generators'] + row['provenance'] + row['references']
            parts.append('<details><summary>Sources, configuration and provenance</summary><p>Physical settings and event counts are recorded in the report and its saved cards. Original run manifests retain their recorded identities.</p><ul>'+''.join('<li>'+link(p,p)+'</li>' for p in dict.fromkeys(details))+'</ul></details></article>')
        parts.append('</section>')
    parts.append('<footer>Regenerate: <code>make report-catalogue</code>. Registry: '+link(REGISTRY,"generated_reports.json")+'. Recorded paths are resolved through the migration registry.</footer><script>document.querySelector("#search").addEventListener("input",e=>{let q=e.target.value.toLowerCase();document.querySelectorAll("article").forEach(a=>a.hidden=!a.textContent.toLowerCase().includes(q));document.querySelectorAll("section").forEach(s=>s.hidden=![...s.querySelectorAll("article")].some(a=>!a.hidden));});</script></main></html>')
    (OUT/"index.html").write_text("\n".join(parts)+"\n")
    (OUT/"catalogue.json").write_text(json.dumps({"schema":1,"reports":entries},indent=2)+"\n")
    for path in links:
        if not path.exists(): raise FileNotFoundError(path)
    return {"collections":len(entries),"documents":sum(len(r['documents']) for r in entries),"checked_links":len(links)}


def write_julia_path_map(registry):
    """Export the shared mapping as stdlib-readable TOML for Julia callers."""
    moves = json.loads((ROOT / "docs/ltspice_paths.json").read_text())["moves"]
    pairs = [(r["old"], r["new"]) for r in moves]
    pairs.extend(("latex_generated/" + slug, "latex_generated/" + spec["group"] + "/" + slug)
                 for slug, spec in registry.items())
    text = "# Generated by make report-catalogue from the two JSON path registries.\n"
    for old, new in pairs:
        text += "\n[[moves]]\nold = " + json.dumps(old) + "\nnew = " + json.dumps(new) + "\n"
    (ROOT / "docs/path_migrations.toml").write_text(text)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--organize", action="store_true")
    parser.add_argument("--clean-intermediates", action="store_true")
    args = parser.parse_args()
    registry = json.loads(REGISTRY.read_text())["reports"]
    write_julia_path_map(registry)
    if args.organize: print("Moved collections:", len(relocate(registry)))
    if args.clean_intermediates: print("Relocated build intermediates:", clean_intermediates(registry))
    print(json.dumps(generate(registry),indent=2))


if __name__ == "__main__":
    main()
