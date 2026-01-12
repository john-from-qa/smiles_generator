# smiles_generator

Python command-line utility that generates molecules in SMILES format from a given core scaffold.

---

## 🔧 Build & run (recommended: Conda) ✅

1. Create and activate the conda environment (uses `conda-forge` RDKit):

```bash
conda create -n smiles python=3.11 -c conda-forge rdkit -y
conda activate smiles
```

2. Run the generator:

```bash
python smiles_generator.py --help
python smiles_generator.py --count 3 --core "C1=CC=C2C(=C1)NC=C2"
```

Output will be written to `random_mol.smi`.

---

## 🧪 Alternative: venv + pip (may fail for RDKit on some platforms)

1. Create and activate a virtual environment:

```bash
python -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

2. Run as above.

---

> ⚠️ Note: If you encounter an error like:
> "A module that was compiled using NumPy 1.x cannot be run in NumPy 2.x",
> pin `numpy<2` in your environment (this repository's `environment.yml` includes that pin).

## 🛠️ Convenience helpers

- `Makefile` targets:
  - `make conda-env` — create/update conda env from `environment.yml`
  - `make pip-env` — create a `.venv` and install `requirements.txt`
  - `make run` — run generator with current Python
  - `make run-conda` — run generator inside the `smiles` conda env

---

## ✅ CI

A GitHub Actions workflow (`.github/workflows/ci.yml`) is included to perform a smoke run on push / PR using the `environment.yml` conda environment.

---

If you'd like, I can also:
- add a `tox` config or a simple test that validates generator output, or
- run a smoke test in this environment now (requires RDKit installed).

Let me know which you'd like me to do next. Thank you!
