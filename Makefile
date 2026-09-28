.PHONY: install test demo lint clean

install:
	python3 -m pip install -e ".[dev,data]"

test:
	python3 -m pytest -q

demo:
	python3 examples/bias_demo.py

clean:
	rm -rf build dist src/*.egg-info .pytest_cache
	find . -name __pycache__ -type d -exec rm -rf {} +
