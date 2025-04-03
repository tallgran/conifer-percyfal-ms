rule save_config:
    """Save configuration report"""
    output:
        report(
            "results/config/{analysis}.config.yaml",
            caption="../report/config.rst",
            category="Configuration",
        ),
    log:
        "logs/results/config/{analysis}.save_config.log",
    script:
        "../scripts/save_config.py"


localrules:
    save_config,
