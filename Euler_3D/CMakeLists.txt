cmake_minimum_required(VERSION 2.8)

project(Euler_3D_2)
enable_language(CXX)
set(CMAKE_CXX_STANDARD 17) # C++17...
set(CMAKE_CXX_STANDARD_REQUIRED ON) #...is required...
set(CMAKE_CXX_EXTENSIONS OFF)

# リリースビルドのためのコンパイラオプション
set(CMAKE_CXX_FLAGS_RELEASE "-DNDEBUG -Wall -Wextra -O3 -flto")

# デバッグビルドのためのコンパイラオプション
set(CMAKE_CXX_FLAGS_DEBUG "-v -Wall -Wextra -g -O0")

set(SRC
        "main.cpp"
        "FluidInitialize.cpp"
        "OutflowBoundary.cpp"
        "ReflectBoundary.cpp"
        )

add_executable(${PROJECT_NAME} ${SRC})
