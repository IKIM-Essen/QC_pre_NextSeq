rule create_sample_sheet:
    output:
        "config/pep/samples.csv",
    params:
        inpath=config["sample-sheet"]["data-path"],
        renaming=config["sample-sheet"]["rename-sample-files"],
    log:
        "logs/create_sample_sheet.log",
    conda:
        "../envs/python.yaml"
    script:
        "../scripts/create_sample_sheet.py"
