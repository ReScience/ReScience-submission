## Introduction


Reference implementation of

Gancarz, G., Grossberg, S. "A Neural Model of the Saccade Generator in the Reticular Formation." 
Neural Networks 11, no. 7-8 (October 1998): 1159-74. doi:10.1016/S0893-6080(98)00096-3.

The model has been implemented in the NEST simulator (2.16.0-beta) using Python (2.7.12) as interface. Extensive documentation can be found on the official website [(http://nest-initiative.org/)](http://nest-initiative.org/).

## NEST Simulator installation instructions

## Linux

1. Create a folder wherein you would like to install NEST (the path to which we refer to by`$NEST`) 

2. Download the latest version of NEST (2.16.0-beta): [https://github.com/nest/nest-simulator/archive/master.zip](https://github.com/nest/nest-simulator/archive/master.zip)

3. Extract the contents of the zip file into the previously created folder (i.e. `$NEST/nest-simulator-master`)

4. Create folders `nest-simulator-master.build` and `nest-simulator-master.install` on the same level as `nest-simulator-master`

5. Configure and build NEST inside the build directory:

```shell
$ cmake -DCMAKE_INSTALL_PREFIX=$NEST/nest-simulator-master.install -DCMAKE-C_COMPILER=gcc 
-DCMAKE_CXX_COMPILER=g++ -Dwith-mpi=ON ../nest-simulator-master
$ make
$ make install
```

6. Add the script `/install/path/bin/nest_vars.sh` to your `.bashrc` file to automatically set paths to the NEST executable and the NEST Python module.

## Windows
Windows and Linux differ considerably. This is the reason why it is difficult to compile NEST natively under Windows. However, it is still possible to use NEST under Windows using one of the following methods:

### Virtual Machines

A virtual machine is a software that lets you use Linux on top of Windows. Popular virtual machine packages are VirtualBox and VMWare. Once you have installed a virtual machine package, you can download OVA images containing NEST and import them into your virtual machine.

### Cygwin

Cywin is a software layer which emulates Unix system calls. NEST should compile and install under Cygwin with the generic installation instructions, but we have not tested this for a long time and do not support NEST under Cygwin at present. Compilation under Cygwin is very slow, but the execution times are comparable to an equivalent Linux installation

## Scripts

Each experiment (simulation) reported in Gancarz & Grossberg (1998) is implemented as a separate Python script.

- Exp1.py
- Exp2.py
- Exp3.py
- Exp4.py
- Exp5.py
- Exp6.py
- Exp7a.py (input as described in original publication)
- Exp7b.py (adjusted input)
- Exp8.py
- Exp9.py


