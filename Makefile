# Simple Makefile to manage environment and run the generator
.PHONY: help conda-env pip-env run run-conda

help:
	@echo "Usage: make <target>"
	@echo "Targets:"
	@echo "  conda-env   Create/update conda environment from environment.yml"
	@echo "  pip-env     Create venv and install pip requirements"
	@echo "  run         Run generator using current Python"
	@echo "  run-conda   Run generator inside the 'smiles' conda env (requires conda)"

conda-env:
	conda env create -f environment.yml -n smiles || conda env update -f environment.yml -n smiles

pip-env:
	python -m venv .venv
	.venv/bin/pip install --upgrade pip
	.venv/bin/pip install -r requirements.txt

run:
	python smiles_generator.py --count 1

run-conda:
	conda run -n smiles python smiles_generator.py --count 1
