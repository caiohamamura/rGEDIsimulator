REM # Download/Install miniconda 
REM # and open this script using Anconda Prompt
call conda create -n gedisimulator -y 
call conda activate gedisimulator

pushd %CONDA_PREFIX%
mkdir src
cd src

call conda install -y git
call conda install -y geotiff hdf5 gdal 
call conda install -y cmake gsl

REM install visual_studio_sdk

git clone https://bitbucket.org/caiohamamura/gedisimulator
git clone https://github.com/caiohamamura/tools
git clone https://github.com/caiohamamura/libclidar

call conda install -y -c menpo wget

wget --no-check-certificate https://www.physics.wisc.edu/~craigm/idl/down/cmpfit-1.2.tar.gz
tar -xf cmpfit-1.2.tar.gz


cd gedisimulator
mkdir build
cd build
cmake -G "Visual Studio 15 2017 Win64" ..
cmake --build . --config Release
copy Release\*.exe %CONDA_PREFIX%\Library\bin\