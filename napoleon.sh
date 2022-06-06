sphinx-apidoc -f -o docs/source .
cd docs
make clean
make html
