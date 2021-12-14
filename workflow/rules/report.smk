rule plot_sample_composition:
    input:
        tsv=rules.generate_krona_input.output.krona_input,
    output:
        chart=report(
            plot_dir / "krona/{sample}.krona.html",
            category="Krona",
            caption=report_dir / "krona.rst",
            ),
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt * GB,
    conda:
        str(env_dir / "krona.yaml")
    log:
        rule_log_dir / "plot_sample_composition/{sample}.log",
    shell:
        "ktImportText {input.tsv} -o {output.chart} &> {log}"


rule composition_report:
    input:
        filter_logs=list(filter_logfiles),
        lineages=list(mykrobe_results)
    output:
        html=report(
            report_dir / "composition.html",
            category="Composition",
            caption=report_dir / "composition.rst",
            ),
    threads: 1
    resources:
        mem_mb=GB,
    conda:
        str(env_dir / "composition_report.yaml")
    params:
        script=scripts_dir / "composition_report.yaml",
        template=report_dir / "composition.html.jinja",
        contam_warning=5.0,
        unmapped_warning=5.0,
    log:
        rule_log_dir / "composition_report.log",
    shell:
        """
        python {params.script} --assignment-dir {params.assignment_dir} \
            --logs-dir {params.logs_dir} \
            --template {params.template} \
            -o {output.html} \
            --contam-warning {params.contam_warning} \
            --unmapped-warning {params.unmapped_warning} 2> {log}
        """


def infer_subsample_log_dirs(logfiles: List[PathLike]) -> List[str]:
    # rule_log_dir / f"subsample_{tech}_reads" / f"{site}" / f"{sample}.log"
    dirs = set()
    for p in logfiles:
        dirs.add(str(Path(p).parent.parent))
    return list(dirs)


rule coverage_report:
    input:
        lineage=expand(
            str(report_dir / "lineage_assignment" / "{sample}.lineage.csv"),
            sample=samples,
        ),
        subsample_logs=list(subsample_logfiles),
    output:
        html=report(
            report_dir / "coverage.html",
            category="Coverage",
            caption=report_dir / "coverage.rst",
            ),
        csv=report(report_dir / "coverage.csv", category="Coverage"),
    threads: 1
    resources:
        mem_mb=GB,
    conda:
        envs["coverage_report"]
    params:
        script=scripts["coverage_report"],
        assignment_dir=lambda wildcards, input: Path(input.lineage[0]).parent,
        logs_dir=lambda wildcards, input: infer_subsample_log_dirs(input.subsample_logs),
    log:
        rule_log_dir / "coverage_report.log",
    shell:
        """
        python {params.script} --assignment-dir {params.assignment_dir} \
            -o {output.html} -c {output.csv} \
            {params.logs_dir} 2> {log}
        """
