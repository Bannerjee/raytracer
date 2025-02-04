cmake_minimum_required(VERSION 3.20)
project(tracer LANGUAGES CXX)
set(CMAKE_CXX_STANDARD 20)


add_executable(
	tracer 
	src/main.cpp
)
target_compile_options(tracer PRIVATE -g -ggdb -pedantic -O3 -Wall -Wextra -Wno-unused-parameter)

