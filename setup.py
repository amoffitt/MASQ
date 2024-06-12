import setuptools


setuptools.setup(
    name="masq",
    url="https://github.com/amoffitt/MASQ",
    packages=setuptools.find_packages(
        where=".", exclude=[
            "tests.*", "tests"
        ],
    ),
    package_data={
        "masq": ["py.typed"],
    },
    scripts=[
    ],
    entry_points="""
    [console_scripts]
    masq_primer_table_to_sd_table=masq.primer_table_to_sd_table:main
    masq_check_loci_plot_and_extend=masq.check_loci_plot_and_extend:main
    masq_combine_within_tag_err=masq.combine_within_tag_err:main
    masq_tag_count_graphs=masq.tag_count_graphs:main
    masq_tag_count_graphs_allregions=masq.tag_count_graphs_allregions:main
    masq_collapse_tags=masq.collapse_tags:main
    masq_combine_reports=masq.combine_reports:main
    masq_extract_variant_info=masq.extract_variant_info:main

    masq_select_enzymes_for_snps=masq.primer_design.select_enzymes_for_snps:main
    """,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
