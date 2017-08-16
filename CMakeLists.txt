project(popkmer)

set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}  -O3  -std=c++0x -Wno-unused-variable -Wno-unused-result -Wno-format-security")

add_executable(popkmer Main.cpp)
