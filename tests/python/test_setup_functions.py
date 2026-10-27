# ---------------------------------
# wisdm - tests for src/setup_functions.py
# ApexRMS, October 2026
# ---------------------------------

import numpy as np
import pandas as pd
import pytest

import setup_functions as sf


# signedDtype ---------------------------------------------------------------

@pytest.mark.parametrize("inputDtype, expectedDtype", [
    (np.uint8, "int16"),
    (np.uint16, "int32"),
    (np.uint32, "int64"),
    (np.int16, "int16"),
    (np.float32, "float32"),
])
def test_signed_dtype_returns_type_that_can_hold_nodata(inputDtype,
                                                       expectedDtype):
    """Unsigned types are widened so the negative nodata value (-9999) fits;
    signed and float types are returned unchanged."""
    assert sf.signedDtype(np.dtype(inputDtype)) == expectedDtype


def test_signed_dtype_rejects_uint64():
    """There is no larger signed integer to promote uint64 to."""
    with pytest.raises(ValueError, match="uint64"):
        sf.signedDtype(np.dtype(np.uint64))


# getNumThreads -------------------------------------------------------------

@pytest.mark.parametrize("enableMultiprocessing, expectedThreads", [
    ("Yes", 4),
    ("No", 1),
])
def test_get_num_threads_follows_multiprocessing_setting(enableMultiprocessing,
                                                         expectedThreads):
    """Uses MaximumJobs only when multiprocessing is turned on."""
    multiprocessingSheet = pd.DataFrame({
        "EnableMultiprocessing": [enableMultiprocessing],
        "MaximumJobs": [4],
    })

    assert sf.getNumThreads(multiprocessingSheet) == expectedThreads
