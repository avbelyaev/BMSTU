распечатать гимпловское представление программы
в gcc-4.9.0 в tree-ssa-dce.c добавляем функцию печати гимпл стейтментов

cd build
make install
//пересобирается gcc

//если make install не сработает то добавляем слудющие временные библиотеки:
//или проипсываем их в bash rc
export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="${CURR}/externals/mpc_install/lib:${LD_LIBRARY_PATH}"
export LD_LIBRARY_PATH="/home/anthony/externals/mpc_install/lib:${LD_LIBRARY_PATH}"

//компилим программу с измененным gcc
cd ../install/bin
./gcc -O2 -o hello hello.c



//еще фиксы для библиотек:
anthony@tokyonyquisd:~/opt_gcc/install/bin$ ./gcc -O2 -o hello hello.c
/home/anthony/opt_gcc/install/libexec/gcc/x86_64-unknown-linux-gnu/4.9.0/cc1: error while loading shared libraries: libmpc.so.2: cannot open shared object file: No such file or directory
anthony@tokyonyquisd:~/opt_gcc/install/bin$ export LD_LIBRARY_PATH=/home/anthony/opt_gcc/externals/mpc_install/lib:${LD_LIBRARY_PATH}
anthony@tokyonyquisd:~/opt_gcc/install/bin$ export LD_LIBRARY_PATH=/usr/lib/x86_64-linux-gnu:${LD_LIBRARY_PATH}
anthony@tokyonyquisd:~/opt_gcc/install/bin$ 
anthony@tokyonyquisd:~/opt_gcc/install/bin$ ./gcc -O2 -o hello hello.c
anthony@tokyonyquisd:~/opt_gcc/install/bin$ 
