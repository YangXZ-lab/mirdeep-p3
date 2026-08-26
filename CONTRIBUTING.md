# Contributing to miRDeep-P3

Thank you for considering contributing to miRDeep-P3! We welcome
bug reports, feature suggestions, documentation improvements, and
pull requests.

## Ways to contribute

- **Report bugs**: open an Issue with the exact command, input files,
  and the full error log
- **Suggest features**: describe the use case and expected behavior
- **Submit pull requests**: code, docs, tests

## Development setup

```bash
git clone https://github.com/YangXZ-lab/mirdeep-p3.git
cd mirdeep-p3

# create environment
conda env create -f mirdp3_environment.yml -n mirdp3
conda activate mirdp3

# optional: download index data for full E2E test
wget https://github.com/YangXZ-lab/mirdeep-p3/releases/download/mirdeep-p3-v3.1.4c-full/data-index.tar.gz
tar xzf data-index.tar.gz -C data/
```

## Before submitting a PR

1. **Run the test suite**:
   ```bash
   ./test.sh
   ```
2. If you change behavior, update or add tests under `tests/`
3. If you change user-facing features, update `README.md` and `docs/`
4. Describe your changes clearly in the PR description

## Code style

- Follow [PEP 8](https://peps.python.org/pep-0008/) for Python code
- Use English for code comments and variable names
- Keep executable scripts with `chmod 755`
- Reference external tools by their canonical names (e.g. `bowtie`, `RNAfold`)

## Licensing

By contributing, you agree that your contributions will be licensed
under the project's BSD 3-Clause License (see `LICENSE`).
