#!/usr/bin/env python3
"""Build config/samples.tsv from the ENA record of BioProject PRJNA472973 (GEO GSE114883).

The study has 24 libraries: 4 donors x {mock-infected, HIV latency model} x {DMSO, SAHA, RMD}.
This workflow analyses the 12 mock-infected libraries, i.e. the drug effect on uninfected
primary CD4+ T cells. Run titles on ENA look like "8_D2_uninf_SAHA".

    python workflow/scripts/make_sample_sheet.py > config/samples.tsv
"""
import csv
import io
import sys
import urllib.request

PROJECT = "PRJNA472973"
FIELDS = "run_accession,sample_alias,sample_title,read_count,fastq_ftp,fastq_md5"
URL = (
    "https://www.ebi.ac.uk/ena/portal/api/filereport"
    f"?accession={PROJECT}&result=read_run&fields={FIELDS}&format=tsv"
)
TREATMENT_ORDER = {"DMSO": 0, "SAHA": 1, "RMD": 2}


def main():
    with urllib.request.urlopen(URL, timeout=120) as response:
        text = response.read().decode()

    rows = []
    for rec in csv.DictReader(io.StringIO(text), delimiter="\t"):
        _, donor, infection, treatment = rec["sample_title"].split("_")
        if infection != "uninf":
            continue
        fq1, fq2 = rec["fastq_ftp"].split(";")
        md5_1, md5_2 = rec["fastq_md5"].split(";")
        rows.append(
            {
                "sample": f"{donor}_{treatment}",
                "donor": donor,
                "treatment": treatment,
                "run": rec["run_accession"],
                "gsm": rec["sample_alias"],
                "read_pairs": rec["read_count"],
                "fastq_1": "https://" + fq1,
                "fastq_2": "https://" + fq2,
                "md5_1": md5_1,
                "md5_2": md5_2,
            }
        )

    if len(rows) != 12:
        sys.exit(f"expected 12 mock-infected libraries, found {len(rows)}")

    rows.sort(key=lambda r: (int(r["donor"][1:]), TREATMENT_ORDER[r["treatment"]]))
    writer = csv.DictWriter(sys.stdout, fieldnames=list(rows[0]), delimiter="\t", lineterminator="\n")
    writer.writeheader()
    writer.writerows(rows)


if __name__ == "__main__":
    main()
