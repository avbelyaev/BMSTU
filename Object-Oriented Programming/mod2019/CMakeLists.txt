cmake_minimum_required(VERSION 3.14)
project(mod2019)

set(CMAKE_CXX_STANDARD 14)

# TASK 1
#add_executable(mod2019 text_analysis/test.cpp text_analysis/textstats.hpp text_analysis/textstats.cpp)
#include_directories(text_analysis)
#configure_file(text_analysis/10.txt 10.txt COPYONLY)

# TASK 2
#add_executable(mod2019 spellcheck/spellcheck.cpp)
#include_directories(spellcheck)
#configure_file(spellcheck/count_big.txt count_big.txt COPYONLY)

# TASK 3
add_executable(mod2019 excel/excel.cpp)
include_directories(excel)