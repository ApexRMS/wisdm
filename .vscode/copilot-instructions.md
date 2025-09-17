# Coding Standards & Guidelines

These instructions guide GitHub Copilot and CodeRabbit in generating and reviewing code for this repository.  
Follow these rules for **style**, **structure**, and **best practices**.

---

## 1. General Principles

- Always write **clear, readable, and maintainable** code.
- Prefer **self-explanatory names** over comments explaining bad names.
- Write code as if the **next maintainer has no context**.
- Avoid unnecessary complexity — **simplicity beats cleverness**.

---

## 2. File & Project Organization

- Use **relative paths** instead of absolute paths.
- Group helper functions in a `helper-functions` file or module.
- Never hard-code paths on local disk; use configuration files or environment variables.
- Organize project folders as:

  - src/ Source code
  - tests/ # Unit and integration tests
  - docs/ # Documentation
  - scripts/ # Utility scripts

- Commit only essential files; ignore build artifacts and large datasets.

---

## 3. Naming Conventions

- **R**: Use `camelCase` for variables and functions (e.g., `processData`).
- **Python**: Use `snake_case` for variables and functions (e.g., `process_data`) unless code is already written using camelCase.
- Name variables as **nouns** and functions as **verbs**.
- Name single objects in singular (e.g., `user`), and collections in plural (e.g., `users`).
- File names: Use kebab-case (lowercase with hyphens) (e.g., `data-cleaning.py`).
- Maintain consistent naming across the project.

---

## 4. Code Formatting

- **Line length**: 80 characters max.
- **Indentation**: 4 spaces for Python, 2 spaces for R.
- Break long function arguments onto multiple lines and align them.
- Separate logical sections of code.
- Avoid trailing whitespace.

---

## 5. Comments & Documentation

- Use docstrings in Python and roxygen2 in R for all functions.
- Comments should explain why, not what.
- Minimize comments as much as possible; object and function names should be self-explanitory

---

## 6. Testing

- **R**: use testthat.
- **Python**: use pytest or unittest.
- Write unit tests for every new function.
- Run tests before committing changes.

---

## 8. Language-Specific Notes

R:

- Assign variables with <- except in function arguments.
- Avoid rgdal and raster — use terra instead.
- Use pipes %>% when possible.

Python

- Follow PEP 8, with the exception of using camelCase
- Use import module as alias (e.g., import pandas as pd), not from module import function.
- Prefer f-strings for formatting.

---

Last updated: 2025-08-21
