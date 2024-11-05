#! /bin/bash -f



## Chech that the PYTHON packages are installed
echo " "
echo " "
echo "INSTALL PYTHON PACKAGES"
echo " "
package_list="numpy agama astropy matplotlib scipy dustmaps pygaia git+https://github.com/keflavich/imf.git"

python3 -m pip install ${package_list}

echo ""
echo "---------------------------------------"
echo "         INSTALLATION DONE"
echo "---------------------------------------"
