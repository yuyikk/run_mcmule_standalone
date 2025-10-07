# Readme
## Compilation && Running
Please modify the following line in the `CMakeLists.txt` and change the path to mcmule installation path
`set(MCMULE_LIB_PATH $ENV{MCMULE_PATH})`

```bash
mkdir build
cd build
cmake ..
make -j
mcmule ./libuser_ep.so < ../card_ep
mcmule ./libuser_ee.so < ../card_ee
mcmule ./libuser_convert_from_fortran.so < ../card_ep
```
## About user files
- `user_ee.cpp` and `user_ep.cpp`: 
these two files are based on the user file in [simple_geant4_package](https://gitlab.com/yannickulrich/simple_geant4_package) but removed all the GEANT4-related codes.
- `user_convert_from_fortran.cpp`:
This file is based on the [fortran version on gitlab](https://gitlab.com/mule-tools/user-library/-/blob/magix/l-p-scattering/prad/prad11_nnlo/user-prad.f95)