"""Resolved paths shared by every Python script in this repository.

Imported by the bootstrap block at the top of each script (`sys.path` already
contains this directory). Override with environment variables:

    GLIOMA_ROOT     repository root                        (default: auto-detected)
    GLIOMA_DATA     raw-data tree                          (default: <sibling>/Datasets)
    GLIOMA_SIBLING  directory holding companion projects   (default: dirname(ROOT))
"""
import os
import pathlib

ROOT = pathlib.Path(__file__).resolve().parent.parent          # repository root

_sibling = os.environ.get("GLIOMA_SIBLING")
SIBLING_ROOT = pathlib.Path(_sibling).resolve() if _sibling else ROOT.parent

_data = os.environ.get("GLIOMA_DATA")
DATA_ROOT = pathlib.Path(_data).resolve() if _data else SIBLING_ROOT / "Datasets"

REF_ROOT = ROOT / "reference"                                   # shipped reference tables

if not (ROOT / "code").is_dir():
    raise SystemExit("repository root does not contain code/ -- set GLIOMA_ROOT to the cloned repository")
