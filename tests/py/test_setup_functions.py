import numpy as np
import pandas as pd
import pytest

from setup_functions import getNumThreads, signedDtype


@pytest.mark.parametrize(
    "dtype, expected",
    [
        (np.dtype("uint8"), "int16"),
        (np.dtype("uint16"), "int32"),
        (np.dtype("uint32"), "int64"),
    ],
)
def test_signedDtype_promotes_unsigned_types(dtype, expected):
    assert signedDtype(dtype) == expected


@pytest.mark.parametrize(
    "dtype", [np.dtype("int16"), np.dtype("int32"), np.dtype("float32"), np.dtype("float64")]
)
def test_signedDtype_leaves_signed_and_float_types_unchanged(dtype):
    assert signedDtype(dtype) == str(dtype)


def test_signedDtype_rejects_uint64():
    with pytest.raises(ValueError, match="uint64"):
        signedDtype(np.dtype("uint64"))


def test_getNumThreads_uses_configured_worker_count_when_multiprocessing_enabled():
    sheet = pd.DataFrame({"EnableMultiprocessing": ["Yes"], "MaximumJobs": [4]})
    assert getNumThreads(sheet) == 4


def test_getNumThreads_defaults_to_one_when_multiprocessing_disabled():
    sheet = pd.DataFrame({"EnableMultiprocessing": ["No"], "MaximumJobs": [8]})
    assert getNumThreads(sheet) == 1
