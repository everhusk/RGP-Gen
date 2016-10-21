RGP Generator
=============

This program was developed to generate custom step size RGP files using the NIST fluid properties.

How to Use
----

RGP File gen is a simple program generating complex files for use with the CFX fluid dynamics package.
Once the program is built (see below about how to build) you may run the program on the command line
and answer the questions with regard to your current project. Be advised that the temperature and pressure
ranges that you enter must cover both the total and static states expected in the whole flow field else
you will end up with range clipping (as will be seen in the CFX log file).

New in this version is the ability to put the answers to the questions in the program arguments. Please be 
aware that the order matters, and the order is the same as the normal question asking order. You will see 
the output however printed in the terminal so you will be made aware of any errors.

How to Build
=====

* Bad news: Refprop is licenced by NIST, and should not be distributed.
* Good news: We have everything sorted for you!

RGP Generator is simple to build, your computer may even have everything required
already installed. Please use the following steps to get setup. Be aware if you have
linux, that you would have to had REFPROP installed on windows to get the code and fluid
files.

1) Copy all fluids from your REFPROP install path (c:/program files/refprop) to the fluids folder
2) Copy all the fortran files to the fortran directory
3) Follow the steps below based on your OS

How to Build (Windows)
-----
1) Install mingw with fortran, if this is forign then watch the video available [here](https://www.youtube.com/watch?v=oVfAU1ziOjg
2) Install cmake from [here](https://cmake.org/download/)
3) Make a folder called build in the current directory
4) Run cmake GUI
* Set `where is the source code located` to the folder with RGP.for 
* Set `where to build the binaries` to the build folder you made 
* Press configure and select `MingGW Makefiles` in `Specify the generator for this project`
* Press generate
5) Open a new command window and ensure the current directory is the build folder
6) Run `mingw32-make` which should be located in the bin directory of where you installed mingw32

How to Build (Linux)
-------

Linux is much nicer!

1) Install fortran (`apt-get install gfortran` or `pacman -S gfortran`) based on your system
2) Install cmake (again using apt or pacman)
3) Navigate to the source folder and make a build directory
4) In the build directory run cmake ..
5) Run make


Licence
=======

[BSD Licence](http://opensource.org/licenses/bsd-license.php)
