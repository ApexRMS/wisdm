# wisdm tests

Unit tests check that small pieces of wisdm code (single functions) still give
the right answer after someone changes the code. They run in seconds and don't
need SyncroSim.

## What's here

```
tests/
├── r/
│   ├── run-tests.R                 # runs every R test
│   ├── setup-source-helpers.R      # loads src/ functions before tests run
│   ├── test-00-helper-functions.R  # tests for src/00-helper-functions.R
│   └── fixtures/                   # small sample data used by tests
└── python/
    ├── conftest.py                 # lets tests import modules from src/
    ├── test_setup_functions.py     # tests for src/setup_functions.py
    └── fixtures/
```

The `tests/` folder sits outside `src/`, so it is not bundled into the
SyncroSim package.

New to testthat or pytest? See [r/FRAMEWORK.md](r/FRAMEWORK.md) and
[python/FRAMEWORK.md](python/FRAMEWORK.md) for a short introduction to each.

## Installing the test packages

R tests need `testthat`. Python tests need `pytest`, `numpy` and `pandas`.

### With conda (recommended)

Install the test packages into the same conda environments wisdm uses, so the
tests run on the same R and Python versions as the package. The environment
names come from `src/wisdm-r-conda.yml` and `src/wisdm-py-conda.yml`; run
`conda env list` to check the names on your machine.

```
conda install -n wisdm-r-conda -c conda-forge --freeze-installed r-testthat
conda install -n wisdm-py-conda -c conda-forge --freeze-installed pytest
```

`numpy` and `pandas` are already in the Python environment.

`--freeze-installed` stops conda from upgrading the packages wisdm has pinned.
If conda reports a conflict, use a separate test environment instead:

```
conda create -n wisdm-test -c conda-forge r-base=4.1.3 r-testthat python=3.12 pytest numpy pandas
```

Then activate the environment before running the tests, e.g.
`conda activate wisdm-r-conda`.

### Without conda

- R: `install.packages("testthat")`
- Python: `pip install pytest numpy pandas`

## Running the tests

Run these from the repository root folder.

**R**

```
Rscript tests/r/run-tests.R
```

**Python**

```
pytest tests/python
```

## When to add a test

We only test what we touch. Add a test when you:

- write a new function
- change what an existing function does
- fix a bug (write a test that would have caught it)

Focus on functions that take inputs and return outputs without reading or
writing SyncroSim datasheets, rasters or files. If a function inside a
transformer script needs testing, move it into a helper functions file first
(e.g. `src/08-fit-model-functions.R`), then test it there.

## How to add a test

1. **Find the test file** for the source file you changed. Test files are named
   after the source file:
   - `src/00-helper-functions.R` → `tests/r/test-00-helper-functions.R`
   - `src/setup_functions.py` → `tests/python/test_setup_functions.py`

   If it doesn't exist yet, create it. R test files must start with `test-`.
   Python test files must start with `test_`. For a new R source file, also
   add a `source()` line for it in `tests/r/setup-source-helpers.R`.

2. **Write the test.** Give it a name that says what should happen, call the
   function with a small input, and check the result.

   R:

   ```r
   test_that("tweak.p replaces 0 and 1 with the nearest in-between values", {
     expect_equal(tweak.p(c(0, 0.2, 0.8, 1)), c(0.2, 0.2, 0.8, 0.8))
   })
   ```

   Python:

   ```python
   def test_get_num_threads_is_one_when_multiprocessing_is_off():
       sheet = pd.DataFrame({"EnableMultiprocessing": ["No"],
                             "MaximumJobs": [4]})
       assert sf.getNumThreads(sheet) == 1
   ```

3. **Work out the expected answer yourself** (by hand or from a trusted
   source), not by copying whatever the function currently returns. Add a
   short comment showing how you got it.

4. **Test the edge cases** that matter: empty input, missing values (`NA` or
   `None`), and inputs that should raise an error.

5. **Run the tests** and make sure they pass before opening a pull request.

## Fixtures

If a test needs more than a few values, put a small CSV in `fixtures/` and
read it in the test. In R, use `test_path("fixtures", "my-file.csv")` to find
it. Keep fixtures tiny (a few rows). Never commit real project data or large
rasters.
