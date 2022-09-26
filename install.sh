if [ ! -f "tables/chspec_CHIANTIv10.fits.gz" ]; then 
    echo -e "Table data not available, downloading it now"
    echo -e "Downloading file from https://www.prl.res.in/ch2xsm/static/chspec_CHIANTIv10.fits.gz"
    echo -e "If you require proxy setting for internet access, ensure that proxy is set for wget in /etc/wgetrc or elsewhere. Otherwise download the file from above link and save it in tables directory and re-run this script."
    wget https://www.prl.res.in/ch2xsm/static/chspec_CHIANTIv10.fits.gz -P ./tables
fi
dir=`pwd`
sed -i 's@#define CHSPEC_TABLEPATH.*@#define CHSPEC_TABLEPATH "'$dir'/tables/"@' readTable.h
ln -sf lmodel_chspec.dat lmodel.dat
echo "initpackage chspec lmodel_chspec.dat `pwd` \nexit" | xspec
echo "lmod chspec . \nexit" | xspec
rm -f *~ *.o *FunctionMap.* lpack_* *.mod Makefile
