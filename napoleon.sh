sphinx-apidoc -f -o docs/source .
cd docs
make clean
make html
rm -f docs/source/setup.rst