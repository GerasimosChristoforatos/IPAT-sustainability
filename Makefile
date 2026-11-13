OWID_COMMIT := 22d9a01321ddec951b36f2bfc66e05dc6ad50ad8
PY := python3

.PHONY: all fetch data test
all: data test

fetch: data/raw/owid-co2-data.csv
data/raw/owid-co2-data.csv:
	curl -sS -L -o $@ https://raw.githubusercontent.com/owid/co2-data/$(OWID_COMMIT)/owid-co2-data.csv

data: fetch
	cd data/raw && sha256sum -c SHA256SUMS
	$(PY) src/load.py

test:
	$(PY) -m unittest discover -s tests -v
