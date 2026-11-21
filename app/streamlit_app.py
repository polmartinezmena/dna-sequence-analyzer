"""Launcher for local use and for Streamlit Community Cloud.

Cloud runs this file from the repository root. The package lives in
``src/``, so that folder is added to the path before importing it.
``dna-analyzer web`` runs this same file when it is available.
"""

import sys
from pathlib import Path

_ROOT = Path(__file__).resolve().parents[1]
_SRC = _ROOT / "src"
if _SRC.is_dir():
    sys.path.insert(0, str(_SRC))

from dna_sequence_analyzer.webapp import main

main()
