rule mykrobe:
    input:
        reads=rules.subsample_reads.output.reads,
    output:
        report=amr_dir / "{sample}.mykrobe.json",
    shadow:
        "shallow"
    resources:
        mem_mb=int(4 * GB),
    container:
        containers["mykrobe"]
    log:
        rule_log_dir / "mykrobe/{sample}.log",
    params:
        opts=" ".join(
            [
                "--force",
                "-A",
                "--debug",
                "--format json",
                "-e 0.08",
                "--min_proportion_expected_depth 0.20",
                "--ploidy haploid",
                "--species tb",
            ]
        ),
    threads: 4
    shell:
        """
        mykrobe predict {params.opts} -o {output.report} -i {input.reads} \
          --sample {wildcards.sample} \
          -t {threads} -m {resources.mem_mb}MB > {log} 2>&1
        """
