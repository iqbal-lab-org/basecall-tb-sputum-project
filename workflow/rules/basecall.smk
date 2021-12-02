def infer_barcode_opts(wildcards, threads: int) -> str:
    run = wildcards.run
    run_df = samplesheet.query("run == @run")
    if len(run_df) == 1:  # singleplex
        return ""
    elif len(run_df) == 0:
        raise KeyError(f"Nanopore run '{run}' was not found in the samplesheet")

    kits = " ".join(set(run_df["barcode_kit"]))
    return f'--barcode_kits "{kits}" --trim_barcodes --num_barcode_threads {threads}'


rule basecall:
    input:
        fast5=fast5_dir / "{run}/",
    output:
        summary=(
            data_dir
            / f"basecalls/guppy_v{GUPPY_VERSION}/{{run}}/sequencing_summary.txt"
        ),
        save_path=directory(data_dir / f"basecalls/guppy_v{GUPPY_VERSION}/{{run}}/"),
    threads: 2
    resources:
        mem_mb=int(4 * GB),
    container:
        containers["guppy-gpu"]
    params:
        extras=" ".join(
            [
                "--recursive",
                "--compress_fastq",
                "--device 'cuda:all:100%'",
                f"-c {config['basecall_config']}",
            ]
        ),
        barcode_opts=infer_barcode_opts,
    log:
        rule_log_dir / f"basecall/guppy_v{GUPPY_VERSION}/{{run}}.log",
    shell:
        """
        guppy_basecaller {params.extras} {params.barcode_opts} \
            --num_callers {threads} \
            -i {input.fast5} \
            -s {output.save_path} > {log} 2>&1
        """
