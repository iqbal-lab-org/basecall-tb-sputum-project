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
        lineages=list(mykrobe_results),
        subsample_logs=list(rasusa_logs),
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
        covg_warning=30,
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


rule amr_report:
    input:
        predictions=expand(amr_dir / "{sample}.mykrobe.json", sample=samplesheet.index),
        template=report_dir / "amr.html.jinja",
    output:
        report=report(
            report_dir / "amr.html",
            category="AMR predictions",
            caption=report_dir / "amr.rst",
        ),
    resources:
        mem_mb=int(4 * GB),
    params:
        script=scripts_dir / "amr_report.py",
    log:
        rule_log_dir / "amr_report.log",
    conda:
        str(env_dir / "amr_report.yaml")
    shell:
        "python {params.script} {input.template} {input.predictions} 2> {log} > {output.report}"
