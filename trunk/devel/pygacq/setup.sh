#!/bin/bash

# prevent PYTHONPATH from grabbing stuff we don't want
unset PYTHONPATH

# create a virtual environment to run GACQ in
source /usr/bin/virtualenvwrapper.sh
mkvirtualenv -p `which python` --no-site-packages gacq_devel

dname=`dirname $0`
gempy=`readlink -f $dname/../../`

# add gemini_python to the PYTHONPATH
add2virtualenv $gempy

# install the other needed bits
pip install numpy
pip install Cython
pip install nose
pip install pyfits

# don't use the fortran compiler in /astro
unset F77
pip install scipy

# scikit for fast image contouring
pip install scikit-image

# matplotlib needs X running, so start DS9, since the tests will need it later
ds9 &
pip install matplotlib

# numdisplay needs to be installed from github fork 
if [ ! -d ${dname}/numdisplay ]
then
  cd ${dname}
  git clone git://github.com/coleb/numdisplay.git
fi
cd numdisplay
git pull
python setup.py install 
cd ..

# download the test data
if [ ! -d ${dname}/test/data ]
then
  git clone git://github.com/coleb/gacqtestdata.git ${dname}/test/data
fi

cd ${dname}
# run the tests without DS9 interaction
GACQUI=fast_test nosetests -v test/

# now run the tests with DS9 display
nosetests -v test/

echo "============================================================="
echo "run 'workon gacq_devel' to enter GACQ development environment"
echo "============================================================="
