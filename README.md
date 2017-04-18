# trajcomp
Trajectory Computing Library

# Installing for R
The trajectory computing library is a set of C++ headers. Additionally, there is an R library.

The R library depedns on some packages and you can prepare your system with

    sudo R --no-save < InstallDependencies.R
    sudo make
   
This will pull all needed packages before compiling and installing the library.


# Using C++ headers
You can just copy the header directory from trajcomp/r-package/trajcomp/src/trajcomp to /usr/include
and use it in all your C++ projects.
