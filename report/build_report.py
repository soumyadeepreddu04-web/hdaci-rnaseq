#!/usr/bin/env python3
"""Compile report/report.typ into report/hdaci_rnaseq_report.pdf.

The Typst source reads its numbers from results/summary.json and the QC tables, so the PDF
always matches the last workflow run. Run from anywhere:

    python report/build_report.py
"""
from pathlib import Path

import typst

ROOT = Path(__file__).resolve().parents[1]

typst.compile(
    str(ROOT / "report" / "report.typ"),
    output=str(ROOT / "report" / "hdaci_rnaseq_report.pdf"),
    root=str(ROOT),
)
print("wrote report/hdaci_rnaseq_report.pdf")
