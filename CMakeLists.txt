CMAKE_MINIMUM_REQUIRED(VERSION 3.20)

PROJECT(MATRIX)
SET(CMAKE_CXX_STANDARD 17)
INCLUDE_DIRECTORIES(./include)
ADD_EXECUTABLE(main src/main.cpp src/csc_smatrix.cpp src/triple_smatrix.cpp)# 
