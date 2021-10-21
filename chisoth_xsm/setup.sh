dir=`pwd`
str="char filename[]  = "
sed -i 's+char filename.*+char filename[]  = "'$dir'/IsotSpec_0p1to30p1KeV_0p01Ebin_0p1MKto20MK_eDens1e10_CHIANTIv10.fits";+g' isothermal.cxx
sed -i 's/\r$//g' isothermal.cxx
ln -sf lmodel_isothermal.dat lmodel.dat
echo "initpackage isoth lmodel_isothermal.dat `pwd` \nexit" | xspec
echo "lmod isoth . \nexit" | xspec

