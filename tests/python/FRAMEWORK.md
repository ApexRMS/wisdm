# How pytest works

[pytest](https://docs.pytest.org/) is the most widely used testing tool for
Python. You write small functions that call your code and check the result
with a plain `assert`. pytest finds them, runs them all and reports which
passed and which failed.

For when and where to add tests in wisdm, see [../README.md](../README.md).

## What happens when you run the tests

`pytest tests/python` does the following:

1. Runs `conftest.py` first (here, it adds `src/` to Python's import path so
   tests can `import setup_functions`).
2. Runs every function whose name starts with `test_` in every file named
   `test_*.py`, then prints a summary.

## A test, piece by piece

```python
import pandas as pd
import setup_functions as sf


def test_get_num_threads_is_one_when_multiprocessing_is_off():
    multiprocessingSheet = pd.DataFrame({"EnableMultiprocessing": ["No"],
                                         "MaximumJobs": [4]})

    numThreads = sf.getNumThreads(multiprocessingSheet)

    assert numThreads == 1
```

- **The function name is the description.** It's shown when the test fails,
  so make it say what *should* happen.
- **Set up, run, check.** Most tests follow this pattern: build a small
  input, call the function, then check the result.
- **`assert`** does the checking. If the condition is false, the test fails
  and pytest shows you the values on each side.

Keep each test function focused on one behaviour. If one `assert` fails, the
rest of that function doesn't run.

## Common checks

| Check | Passes when |
|---|---|
| `assert x == y` | `x` equals `y` |
| `assert x` / `assert not x` | `x` is true / false |
| `assert x is None` | `x` is `None` |
| `assert x in y` | `x` is an item in (or a substring of) `y` |
| `assert len(x) == 3` | `x` has 3 items |
| `assert isinstance(x, pd.DataFrame)` | `x` is a DataFrame |
| `assert x == pytest.approx(0.604, abs=1e-3)` | decimal `x` is within 0.001 of 0.604 |
| `with pytest.raises(ValueError):` | the code inside raises a `ValueError` |
| `with pytest.warns(UserWarning):` | the code inside gives a warning |
| `np.testing.assert_allclose(x, y)` | numpy arrays `x` and `y` match |
| `pd.testing.assert_frame_equal(x, y)` | pandas tables `x` and `y` match |

Always use `pytest.approx` for decimals: computers store them with tiny
rounding errors, so `0.1 + 0.2 == 0.3` is `False` in Python. For arrays and
tables, `==` checks item by item, so use the numpy and pandas helpers above.

### Checking for errors

Put the code that should fail inside a `with pytest.raises(...)` block. The
optional `match` checks the error message:

```python
def test_signed_dtype_rejects_uint64():
    with pytest.raises(ValueError, match="uint64"):
        sf.signedDtype(np.dtype(np.uint64))
```

## Running the same test with different inputs

`@pytest.mark.parametrize` runs one test once per set of inputs. Each set is
reported as its own test, so you can see exactly which input failed:

```python
@pytest.mark.parametrize("inputDtype, expectedDtype", [
    (np.uint8, "int16"),
    (np.float32, "float32"),
])
def test_signed_dtype_returns_type_that_can_hold_nodata(inputDtype,
                                                       expectedDtype):
    assert sf.signedDtype(np.dtype(inputDtype)) == expectedDtype
```

## pytest fixtures

In pytest, a "fixture" is a function that prepares something a test needs.
(This is different from the `fixtures/` folder, which holds sample data
files.) Mark the function with `@pytest.fixture`, then list its name as a
test argument and pytest passes in the result:

```python
@pytest.fixture
def multiprocessingSheet():
    return pd.DataFrame({"EnableMultiprocessing": ["Yes"],
                         "MaximumJobs": [4]})


def test_get_num_threads_uses_maximum_jobs(multiprocessingSheet):
    assert sf.getNumThreads(multiprocessingSheet) == 4
```

Fixtures shared by several test files go in `conftest.py`. pytest also has
built-in fixtures. The most useful is `tmp_path`: add it as a test argument
to get a fresh temporary folder that is cleaned up automatically, e.g. for
code that writes files.

To read a file from the `fixtures/` folder, build the path from the test
file's location:

```python
fixturesDir = os.path.join(os.path.dirname(__file__), "fixtures")
siteData = pd.read_csv(os.path.join(fixturesDir, "site-data.csv"))
```

## Running tests

From the repository root:

| To run | Command |
|---|---|
| All Python tests | `pytest tests/python` |
| One file | `pytest tests/python/test_setup_functions.py` |
| One test | `pytest tests/python/test_setup_functions.py::test_signed_dtype_rejects_uint64` |
| Tests whose name contains a word | `pytest tests/python -k num_threads` |

Add `-v` to list each test as it runs, or `-x` to stop at the first failure.

Passing tests show as dots. When a check fails, pytest shows the failing line
and the actual values:

```
    def test_get_num_threads_is_one_when_multiprocessing_is_off():
        ...
>       assert numThreads == 1
E       assert 4 == 1
```

## Going further

Features not covered here, for when you need them:

- [How to invoke pytest](https://docs.pytest.org/en/stable/how-to/usage.html):
  all the command-line options
- [Fixtures](https://docs.pytest.org/en/stable/how-to/fixtures.html): scopes,
  clean-up, sharing fixtures through `conftest.py`
- [Parametrizing tests](https://docs.pytest.org/en/stable/how-to/parametrize.html)
- [Temporary files and folders](https://docs.pytest.org/en/stable/how-to/tmp_path.html):
  using `tmp_path`
- [Monkeypatching and mocking](https://docs.pytest.org/en/stable/how-to/monkeypatch.html):
  replace a function (e.g. a pysyncrosim call) or environment variable
  during a test
- [Skip and xfail](https://docs.pytest.org/en/stable/how-to/skipping.html):
  tests that can't run on some machines, or known failures
- [Plugins](https://docs.pytest.org/en/stable/how-to/plugins.html), e.g.
  [pytest-cov](https://pytest-cov.readthedocs.io/) for test coverage reports
- [API reference](https://docs.pytest.org/en/stable/reference/reference.html):
  every function, including `pytest.approx` and `pytest.raises`
