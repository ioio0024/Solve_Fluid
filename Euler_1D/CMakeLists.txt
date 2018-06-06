cmake_minimum_required(VERSION 2.8)

project(Euler_1D)

set(CMAKE_CXX_STANDARD 17) # C++17...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF)

set(SRC "main.cpp"
        "CipCsl3.cpp"
        "Fluid1D.cpp"
        "FluidInitialize.cpp"
        "Advection.cpp"
    )

add_executable(${PROJECT_NAME} ${SRC})
