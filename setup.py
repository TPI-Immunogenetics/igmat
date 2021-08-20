from setuptools import setup

# Metadata goes in setup.cfg. These are here for GitHub's dependency graph.
setup(
    name="IgMAT",
    install_requires=[
        'prettytable',
        'biopython',
        'pyyaml'
    ],
    extras_require={
    },
)