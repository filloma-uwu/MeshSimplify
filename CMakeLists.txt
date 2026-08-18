cmake_minimum_required(VERSION 3.24)

project(MeshSimplify VERSION 0.1.0 LANGUAGES CXX)

option(MESH_SIMPLIFY_BUILD_TESTS "Build MeshSimplify tests" ON)

add_library(mesh_simplify STATIC
    src/primitive_mesh_analyzer.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.engine.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.offset.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.rectclip.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.triangulation.cpp
)

target_include_directories(mesh_simplify
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/clipper2/CPP/Clipper2Lib/include
)

target_compile_features(mesh_simplify PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(mesh_simplify PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(mesh_simplify PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(mesh-simplify-analyze tools/primitive_mesh_analyze.cpp)
target_link_libraries(mesh-simplify-analyze PRIVATE mesh_simplify)
target_compile_features(mesh-simplify-analyze PRIVATE cxx_std_20)

if(MSVC)
    target_compile_options(mesh-simplify-analyze PRIVATE /W4 /permissive- /utf-8)
else()
    target_compile_options(mesh-simplify-analyze PRIVATE -Wall -Wextra -Wpedantic)
endif()

if(MESH_SIMPLIFY_BUILD_TESTS)
    enable_testing()

    add_executable(mesh_simplify_analyzer_test
        tests/primitive_mesh_analyzer_test.cpp
    )
    target_link_libraries(mesh_simplify_analyzer_test PRIVATE mesh_simplify)
    target_compile_features(mesh_simplify_analyzer_test PRIVATE cxx_std_20)
    add_test(NAME mesh_simplify_analyzer_test COMMAND mesh_simplify_analyzer_test)
endif()
