#wget https://www.prl.res.in/ch2xsm/static/chspec_CHIANTIv10.fits -P ./tables
dir=`pwd`
sed -i 's@#define CHSPEC_TABLEPATH.*@#define CHSPEC_TABLEPATH "'$dir'/tables/"@' readTable.h
ln -sf lmodel_chspec.dat lmodel.dat
echo "initpackage chspec lmodel_chspec.dat `pwd` \nexit" | xspec
echo "lmod chspec . \nexit" | xspec
rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile
