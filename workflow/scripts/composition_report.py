"""This script generates a HTML file with a table containing information about the
composition and lineage of each sample.
"""
import sys

sys.stderr = open(snakemake.log[0], "w")
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import TextIO, Tuple
import json
import jinja2
import pandas as pd

NORD0 = "#2e3440"
NORD11 = "#bf616a"
NORD12 = "#d08770"


class RipgrepError(Exception):
    pass


def ripgrep_search(file: Path) -> Tuple[int, int, int]:
    pattern = "\[INFO\]:\s(?P<num>\d+)\sread"
    extra_params = ["--replace", "$num", "--only-matching", "--no-line-number"]
    process = subprocess.Popen(
        ["rg", *extra_params, pattern, str(file)],
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        encoding="utf-8",
    )
    exit_code = process.wait()
    if exit_code != 0:
        raise RipgrepError(
            f"Failed to execute the rg command on the file {file} due to the "
            f"following error:\n{process.stderr.read()}"
        )
    values = map(int, map(str.rstrip, process.stdout.readlines()))
    to_keep, contam, unmapped = list(values)
    return to_keep, contam, unmapped


contam_warning = snakemake.params.get("contam_warning", 5.0) / 100
unmapped_warning = snakemake.params.get("unmapped_warning", 5.0) / 100


def highlight_high_contam(col: pd.Series):
    """Highlights cells in a column if their contamination level is over the
    threshold
    """
    return [
        f"background-color: {NORD12}" if val > contam_warning else "" for val in col
    ]


def highlight_high_unmapped(col: pd.Series):
    """Highlights cells if their unmapped level is over the threshold"""
    return [
        f"background-color: {NORD12}" if val > unmapped_warning else "" for val in col
    ]


def highlight_abnormal_lineages(col: pd.Series):
    """Highlights cells if their lineage is not one of the numbered majors."""
    return [f"background-color: {NORD12}" if v.isalpha() else "" for v in col]


data = defaultdict(dict)
logfiles = snakemake.input.filter_logfiles
for file in logfiles:
    sample = file.name.split(".")[0]
    num_keep, num_contam, num_unmapped = ripgrep_search(file)
    total = sum([num_keep, num_contam, num_unmapped])
    data[sample].update(
        {
            f"keep": num_keep,
            f"keep%": num_keep / total,
            f"contam": num_contam,
            f"contam%": num_contam / total,
            f"unmapped": num_unmapped,
            f"unmapped%": num_unmapped / total,
            f"total": total,
        }
    )

assignment_files = snakemake.input.lineages
for file in assignment_files:
    sample = file.name.split(".")[0]
    with open(file) as fp:
        lineage = json.load(fp)[sample]["phylogenetics"]["lineage"]["lineage"][0]
        lineage = lineage.replace("lineage", "")
        data[sample]["lineage"] = lineage

df = pd.DataFrame(data).T
df.index.name = "sample"

percent_format_cols = [s for s in data[list(data.keys())[0]].keys() if s.endswith("%")]
df_styled = (
    df.style.apply(highlight_high_contam, subset=["contam%"])
    .apply(highlight_abnormal_lineages, subset=["lineage"])
    .apply(highlight_high_unmapped, subset=["unmapped%"])
    .format("{:.2%}", subset=percent_format_cols)
)

table_html = df_styled.render()
template = snakemake.params.template
template_content = Path(template).read_text()
html = jinja2.Template(template_content).render(
    table=table_html,
    contam_warning=contam_warning,
    unmapped_warning=unmapped_warning,
)
outfile = snakemake.output.html
outfile.write(html)
outfile.close()
