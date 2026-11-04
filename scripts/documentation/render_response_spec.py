#!/usr/bin/env python3
"""Render the joint response production specification for browser review."""
from pathlib import Path
from markdown_it import MarkdownIt

ROOT = Path(__file__).resolve().parents[2]
source = ROOT / "docs/response_production_spec.md"
target = source.with_suffix(".html")
body = MarkdownIt("commonmark").enable("table").render(source.read_text())
style = """body{font:16px/1.55 system-ui,sans-serif;color:#243348;max-width:1100px;
margin:32px auto;padding:0 22px}h1,h2{line-height:1.2}h1{margin-bottom:1.4em}
h2{border-bottom:1px solid #dce3ec;padding-bottom:.25em;margin-top:2em}
table{border-collapse:collapse;width:100%;font-size:14px}th,td{border:1px solid #dce3ec;
padding:8px;text-align:left;vertical-align:top}th{background:#edf2f6}
code{background:#f3f6f8;padding:1px 3px;overflow-wrap:anywhere}
ol li{margin-bottom:.45em}"""
target.write_text('<!doctype html><html lang="en"><meta charset="utf-8">'
                  '<title>Ion and SiPM response production specification</title>'
                  f'<style>{style}</style>{body}</html>')
print(target)
