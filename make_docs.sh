
echo ''
echo '         PIP INSTALL MODULE         '
echo '------------------------------------'
pip install .
cd docs
echo ''
echo '            MAKE  CLEAN              '
echo '------------------------------------'
make clean
echo ''
echo '            MAKE  HTML              '
echo '------------------------------------'
make html 
