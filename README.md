# Project Filius #
----

## Fixed Income and Liability Instruments Utility Set ##

### Summary ###

* Core feature  
	Fixed income derivatives pricing framework and market calibration tool

* Description  
    * API assessment style is designed to be single threaded access with configuration based interation paradigm
    * Performance is ensured by hyper parameters that manage multithreaded calculations internally  
    * If you want to integrate this model to your multithreaded system, either/or:  
        (1) Do not share model object between threads, or use separate model objects for each thread   
        (2) Lock model API calls to ensure sequential assess

* Currently supported interest rate models  
    1. Black  
    2. Orstein-Ulhenback system with multiple degrees of freedom  
    * Can easily be extended to integrate support for additional interest rate models  

* Currently supported data calibration method  
    1. Simulated Annealing with Gradient Descent (Global Optimizer)  

* Future support/To do:  
    * Feed-Forward Convolutional Neural Network Pricing Support  
    * Statistical Model Selection Module

* Theoretical references  
    [G2++](doc/G2++.pdf)  
    [Deep Learning](http://www.deeplearningbook.org/)  

* Ideas and designs by author  
    [Model Design - Work in progress](doc/G2++_Math.pdf)  

* Version  
	1.0

### How do I get set up? ###

* Dependencies  
    C++17
* How to build project  
    make (compiles program)  
    make build (creates ojb directories per machine os and architecture, then compiles program, do this if you are compling project for the first time)  
    make clean (clean obj and lib directories for clean build)
* How to run tests  
    make run
* Supported additional test utilities  
    gdb  
    valgrind

### Compiler options ###
* \_\_DEBUG\_\_: Pass this macro flag to compiler will lead to more verbose debug information whenever it is available
* \_\_REGEN\_\_: Pass this macro flag to compiler will force recalculation with refreshed random numbers in in-bulk getZCBP method

### Who do I talk to? ###

* Repo owner  
	Author: James Ding  
    Email: james.ding.illinois@gmail.com

### Usage in academic research or non-profit project ###

This project is developed purely on voluntary basis.
To support the author, please acknowledge this work when you use it in your project.
You can do so by citing this github repo.  
And also feel free to email repo owner if you would like to comment on or contribute to this project.

### Instruction ###

Project Filius is a free software project licensed under [GNU General Public License 3](LICENSE).  
Filius  Copyright (C) 2016-2018  James Ding  
It is provided as is WITHOUT ANY WARRANTY.
You can redistribute it and/or modify it under the terms of the GNU General Public License 3 published by the Free Software Foundation.
