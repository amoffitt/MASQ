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

    """,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Programming Language :: Python :: 3.9",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.9",
)
