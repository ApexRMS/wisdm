# wisdm test suite

## Layout

```
tests/
  R/
    testthat.R                              # entrypoint: Rscript tests/R/testthat.R
    testthat/
      helper-syncrosim.R                    # shared helpers (not a test file itself)
      test-helper-functions.R               # unit tests, src/00-helper-functions.R
      test-data-prep-functions.R            # unit tests, src/05-training-testing-data-prep-functions.R
      test-integration-template-library.R   # builds the "WISDM Example" online template, runs it
      test-integration-synthetic-library.R  # builds a library from scratch, runs it
  py/
    conftest.py
    test_setup_functions.py                 # unit tests, src/setup_functions.py
```

## Prerequisites

- **Unit tests** (`test-helper-functions.R`, `test-data-prep-functions.R`, `test_setup_functions.py`)
  need only a plain R/Python install with the packages the tested functions themselves use
  (base R + `methods`; Python + `numpy`/`pandas`). No SyncroSim installation required.
- **Integration tests** (`test-integration-*.R`) additionally require SyncroSim Studio, the
  `wisdm` package installed from the SyncroSim server, and its conda environments
  (`src/wisdm-r-conda.yml`, `src/wisdm-py-conda.yml`) set up, since they drive the real
  transformer pipeline through `rsyncrosim`. They self-skip (via `skip_if_not(syncrosim_available())`)
  when no SyncroSim session can be reached, so it's safe to run the whole suite on a machine
  without SyncroSim installed.

R package requirements: `testthat`, `rsyncrosim`, `terra` (for the synthetic-library test).

## Running

From the repo root:

```sh
# R suite
Rscript tests/R/testthat.R

# Python suite
pytest tests/py -v
```

## Runtime expectations

- Unit tests: seconds.
- Integration tests: roughly 1-5 minutes each. Both restrict the pipeline (via the
  `core_Pipeline` datasheet, see `wisdm_fast_pipeline_stages` in `helper-syncrosim.R`) to data
  prep + a single GLM fit + spatial prediction, skipping hyperparameter tuning, the other 5
  model types, and the ensemble step -- this is intentional and keeps the suite fast rather
  than being full end-to-end coverage of every model type.

## Notes for maintainers

- `test-integration-template-library.R` depends on the wisdm package's hosted **"WISDM Example"**
  online template (see `docs/getting_started.md`) staying available under that exact name, and on
  its scenario still being named **"Brewer's Sparrow"**. If either is renamed, update the test.
- `test-integration-synthetic-library.R` has no such dependency -- it's the fallback to use when
  validating changes offline or before a new template/scenario name has propagated to the
  SyncroSim server.
- Both integration tests were authored without a local SyncroSim install available to run them
  against, based on the `rsyncrosim` API and `src/package.xml` datasheet/transformer definitions.
  Confirm `ssimLibrary(..., template = ...)` and the `core_Pipeline` datasheet columns
  (`StageNameId`, `RunOrder`) against the actual installed `rsyncrosim` version the first time
  these run, and adjust if the API has moved on.
