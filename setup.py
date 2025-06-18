from setuptools import setup, find_packages

setup(
    name="read-break",
    version="0.1.0",
    packages=find_packages(where = "src"),
    package_dir={"": "src"},
    python_requires=">=3.9",
    install_requires=[],
entry_points={
            "console_scripts": [
                "read_break=main:main",
            ],
    },
    author="Dan Levy",
    author_email="levy@cshl.edu",
    description="Structured pipeline for parsing sequencing reads using Jinja",
)
