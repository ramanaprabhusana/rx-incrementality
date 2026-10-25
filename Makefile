.PHONY: install test demo data data-anticoagulant reproduce estimate sensitivity figures replicate provenance clean

install:
	python3 -m pip install -e ".[dev]"

test:
	python3 -m pytest -q

demo:
	python3 examples/bias_demo.py

# Downloads public CMS data into data/raw (about 1.5 GB, roughly 45 minutes).
data:
	python3 scripts/fetch_class_panel.py --class diabetes --start 2019 --end 2024
	./scripts/stream_all_years.sh 2019 2020 2021 2022 2023 2024

# The anticoagulant replication's inputs (about 1 GB more).
data-anticoagulant:
	python3 scripts/fetch_class_panel.py --class anticoagulant --start 2019 --end 2024
	CLASS=anticoagulant ./scripts/stream_all_years.sh 2019 2020 2021 2022 2023 2024

# Everything downstream of the data, in dependency order.
reproduce: estimate sensitivity figures replicate provenance

estimate:
	python3 scripts/estimate_real.py

sensitivity:
	python3 scripts/sensitivity.py

figures:
	python3 scripts/make_figures.py

replicate:
	python3 scripts/replicate_class.py --class anticoagulant

provenance:
	python3 scripts/provenance.py

clean:
	rm -rf build dist src/*.egg-info .pytest_cache
	find . -name __pycache__ -type d -exec rm -rf {} +
