"""Shim so that `pip install -e .` works with pip older than 21.3.

Everything is declared in pyproject.toml; this file exists only because older
pip cannot do an editable install from PEP 621 metadata alone. Installing is
optional: the Makefile puts src on PYTHONPATH.
"""
from setuptools import setup

setup()
