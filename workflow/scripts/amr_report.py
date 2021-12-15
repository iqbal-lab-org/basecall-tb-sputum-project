import json
import sys
from collections import defaultdict
from pathlib import Path

import jinja2


def extract_mykrobe_predictions(fp):
    data = json.load(fp)
    sample = list(data.keys())[0]
    return data[sample]["susceptibility"]


def evidence2html(ev: dict) -> str:
    html = []
    for variant in ev:
        info = ev[variant]["info"]
        ex_dp = info["expected_depths"][0]
        ref_cov = info["coverage"]["reference"]["median_depth"]
        alt_cov = info["coverage"]["alternate"]["median_depth"]
        s = f"Variant: {variant}<br>Expected|Ref|Alt coverage: {ex_dp}|{ref_cov}|{alt_cov}"
        html.append(s)
    return "<hr>".join(html)


def main():
    files = sys.argv[2:]
    rows = []
    evidence = defaultdict(dict)
    for f in map(Path, files):
        sample = f.name.split(".")[0]
        with open(f) as fp:
            predictions = extract_mykrobe_predictions(fp)
        row = [sample]
        row.extend(
            [predictions[drug]["predict"] for drug in sorted(predictions.keys())]
        )
        rows.append(row)
        for drug in predictions:
            if (pred := predictions[drug])["predict"] == "R":
                ev = pred["called_by"]
                evidence[sample][drug] = evidence2html(ev)

    header = ["Sample"]
    header.extend(sorted(predictions.keys()))

    template = sys.argv[1]
    template_content = Path(template).read_text()
    html = jinja2.Template(template_content).render(
        header=header, rows=sorted(rows), evidence=evidence
    )
    outfile = sys.stdout
    outfile.write(html)


if __name__ == "__main__":
    main()
