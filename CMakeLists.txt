cmake_minimum_required(VERSION 3.20)
project(IDrive LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)


add_executable(
	demo 
	src/main.cpp
)
target_compile_options(demo PRIVATE -g -ggdb -pedantic -O0 -Wall -Wextra -Wno-unused-parameter)

