build_env:
	conda remove --yes -n scins_env --all;
	conda create --yes -n scins_env -c conda-forge python==3.11.8
	conda run -n scins_env pip install -r requirements.txt
	conda run -n scins_env pip install -i https://pypi.python.org/simple/ .

clean:
	rm -fr build/
	rm -fr scins.egg-info
	rm -fr scins.egg

install: clean
	pip install .

test:
	conda run -n scins_env python -m unittest tests/test_scins.py
