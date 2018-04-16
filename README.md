## Introduction


Reference implementation of

Gancarz, G., Grossberg, S. "A Neural Model of the Saccade Generator in the Reticular Formation." 
Neural Networks 11, no. 7-8 (October 1998): 1159-74. <a href="https://doi.org/10.1016/S0893-6080(98)00096-3">doi:10.1016/S0893-6080(98)00096-3</a>.

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

6. Source the script `$NEST/nest-simulator-master.install/bin/nest_vars.sh` in your `.bashrc` file to automatically set paths to the NEST executable and the NEST Python module.

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

## Platform Information

### OS & Architecture
Linux #43~16.04.1-Ubuntu SMP Wed Mar 14 17:48:43 UTC 2018 x86_64 GNU/Linux
### Hardware 
flags: fpu vme de pse tsc msr pae mce cx8 apic sep mtrr pge mca cmov pat pse36 clflush dts acpi mmx fxsr sse sse2 ss ht tm pbe syscall nx pdpe1gb rdtscp lm constant_tsc arch_perfmon pebs bts rep_good nopl xtopology nonstop_tsc cpuid aperfmperf pni pclmulqdq dtes64 monitor ds_cpl vmx smx est tm2 ssse3 sdbg fma cx16 xtpr pdcm pcid dca sse4_1 sse4_2 x2apic movbe popcnt tsc_deadline_timer aes xsave avx f16c rdrand lahf_lm abm 3dnowprefetch cpuid_fault epb cat_l3 cdp_l3 invpcid_single pti retpoline intel_ppin intel_pt tpr_shadow vnmi flexpriority ept vpid fsgsbase tsc_adjust bmi1 hle avx2 smep bmi2 erms invpcid rtm cqm rdt_a rdseed adx smap xsaveopt cqm_llc cqm_occup_llc cqm_mbm_total cqm_mbm_local dtherm ida arat pln pts
model name: Intel(R) Xeon(R) CPU E5-1650 v4 @ 3.60GHz
vendor_id: GenuineIntel
### C Compiler
g++ (Ubuntu 5.4.0-6ubuntu1~16.04.9) 5.4.0 20160609
Copyright (C) 2015 Free Software Foundation, Inc.
This is free software; see the source for copying conditions.  There is NO
warranty; not even for MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
### Cmake
version 3.5.1
### NEST
version 2.16.0-beta
### Python 2.7
Python: 2.7.14 |Anaconda, Inc.| (default, Oct 16 2017, 17:29:19) 
[GCC 7.2.0]

NumPy: 1.13.3

SciPy: 0.19.1

matplotlib: 2.1.0
### Python 3.5
Python: 3.5.2 (default, Nov 23 2017, 16:37:01) 
[GCC 5.4.0 20160609]

NumPy: 1.11.0

SciPy: 0.17.0

matplotlib: 1.5.1

