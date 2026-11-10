# ---------------------------------
# wisdm - test setup
# ApexRMS, October 2026
# ---------------------------------

# pytest runs this file before the tests. It adds src/ to the import path so
# tests can import wisdm modules directly (e.g. `import setup_functions`).

import os
import sys

srcDir = os.path.abspath(
    os.path.join(os.path.dirname(__file__), "..", "..", "src"))
sys.path.insert(0, srcDir)
