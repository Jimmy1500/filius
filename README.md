# Project Filius #
----

## Fix Income and Liability Instruments Utility Set ##

### Summary ###

* Core feature  
	Fix income derivative pricing framework and market calibration tool

* Description  
    Currently supports Gaussian Orstein-Ulhenback interest rate model with multiple degrees of freedom (G2++)  
    Can be easily be extended to integrate support for additional interest rate models  
    Performance is ensured by hyper parameters that manages multithreaded calculations internally 

* Version  
	0.1

### How do I get set up? ###

* Dependencies  
    C++17
* How to run tests  
    gdb  
    valgrind

### Who do I talk to? ###

* Repo owner  
	Owner: Jimmy1500  
    Email: james.ding.illinois@gmail.com

### Usage in academic research or non-profit project ###

This project is developed completely on voluntary basis.
To support our developers, please acknowledge our work when you use Filius in your project.
You can do so by citing this github repo.  
And also feel free to email repo owner if you would like to comment on or contribute to this project.

### Instruction ###

Project Filius is a free software project licensed under [GNU General Public License 3](LICENSE).
It is provided as is WITHOUT ANY WARRANTY.
User will be liable for any consequence from using, redistributing or modifying this software.

### Compiler options ###
* \_\_DEBUG\_\_: Pass this macro flag to compiler will lead to more verbose debug information whenever it is available
* \_\_REGEN\_\_: Pass this macro flag to compiler will force recalculation with refreshed random numbers in in-bulk getZCBP method
