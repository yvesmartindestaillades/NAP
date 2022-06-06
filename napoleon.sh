echo ''
echo '            NAPOLEON                '
echo '------------------------------------'
sphinx-apidoc -f -o docs/source .
cd docs
echo ''
echo '            MAKE  CLEAN              '
echo '------------------------------------'
make clean
echo ''
echo '            MAKE  HTML              '
echo '------------------------------------'
make html 
rm -f docs/source/setup.rst