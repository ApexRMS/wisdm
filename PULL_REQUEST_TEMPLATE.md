## What does this change?

<!-- 1-3 sentences in plain language. Link the issue if there is one,
e.g. "Fixes #12". -->

## How did you check that it works?

<!-- e.g. "Ran the R tests", "Ran the example library in SyncroSim and the
maps look the same as before". -->

## Checklist

<!-- Tick each box that is true. If an item doesn't apply to your change,
write "N/A" after it. Reviewers are happy to help with anything you're
unsure about. -->

### Tests

See [tests/README.md](https://github.com/ApexRMS/wisdm/blob/main/tests/README.md)
for when and how to add tests.

- [ ] I added or updated a unit test for each function I added or changed.
- [ ] If this fixes a bug, I added a test that would have caught it.
- [ ] My tests check answers I worked out myself, not just whatever the code
      currently returns.
- [ ] New logic is in a helper functions file (not inside a transformer
      script), so it can be tested.
- [ ] All tests pass: `Rscript tests/r/run-tests.R` and/or `pytest tests/python`.
- [ ] Any new files in `fixtures/` are small and contain no real project data.

### Easy to read

- [ ] Each new script starts with a short header: name, organization, date
      and what the script does.
- [ ] Each new function has a short description (roxygen comment in R,
      docstring in Python).
- [ ] Names say what things are: variables are nouns, functions are verbs, and
      they match the style of the file I edited (wisdm uses camelCase).
- [ ] Comments explain *why* the code does something, not *what* it does.
- [ ] Code that was repeated is turned into a function.

### Organized

- [ ] Packages, input files and settings are loaded at the top of the script.
- [ ] The script is split into labelled sections.

### Works on other computers

- [ ] No file paths that only exist on my computer (use relative paths).
- [ ] No numbers or text typed directly into the code that a user might need
      to change. Put them at the top of the script, or in `src/00-constants.R`
      or `src/setup_functions.py` if they're shared.
- [ ] New output files are named like the existing ones (e.g.
      `example-filename.tif`).
