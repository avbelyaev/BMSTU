
- установить Ubuntu 14.04 LTS 64-bit
- развернуть архив gcc_how_to_build.zip
- загрузить gcc-4.9.0.tar.bz2 (gcc.gnu.org), положить его рядом c ./_build_it.sh и развернуть
- выполнить ./_apt_get.sh	
- войти в ./externals
- развернуть все архивы
- выполнить
	./_build_gmp.sh
	./_build_ppl.sh
	./_build_mpfr.sh
	./_build_mpc.sh
	./_build_cloog.sh
- cd ..
- выполнить ./_build_it.sh
- в ./install должен лежать готовый компилятор


