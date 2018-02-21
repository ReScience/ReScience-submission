```shell
$ cmake -DCMAKE_INSTALL_PREFIX=$NEST/nest-simulator-master.install -DCMAKE-C_COMPILER=gcc 
-DCMAKE_CXX_COMPILER=g++ -Dwith-mpi=ON ../nest-simulator-master
$ make
$ make install
```

`$NEST` refers to the path to the current folder (`../nest`).
