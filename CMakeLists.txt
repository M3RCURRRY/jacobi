cmake_minimum_required(VERSION 3.17)
project(jcb LANGUAGES C)

find_package(OpenSSL REQUIRED)
add_executable(jcb jacobi.c jacobi.h main.c)
target_link_libraries(jcb OpenSSL::Crypto)