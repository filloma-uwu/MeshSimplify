cmake_minimum_required(VERSION 3.24)

project(PQSSProxyMesh VERSION 0.1.0 LANGUAGES CXX)

option(PQSS_PROXY_MESH_BUILD_TESTS "Build PQSSProxyMesh tests" ON)
option(PQSS_PROXY_MESH_BUILD_REAL_SCENE_COMPARE "Build the PQSS real-scene comparison tool" ON)
set(PQSS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../PQSS" CACHE PATH "Path to the PQSS source tree")

add_library(pqss_proxy_mesh STATIC
    src/mesh_pool.cpp
    src/options.cpp
    src/primitive_mesh_analyzer.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.engine.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.offset.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.rectclip.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.triangulation.cpp
)

target_include_directories(pqss_proxy_mesh
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/clipper2/CPP/Clipper2Lib/include
)

target_compile_features(pqss_proxy_mesh PUBLIC cxx_std_20)

if(MSVC)
    target_compile_options(pqss_proxy_mesh PRIVATE /W4 /permissive-)
else()
    target_compile_options(pqss_proxy_mesh PRIVATE -Wall -Wextra -Wpedantic)
endif()

add_executable(pqss-proxy-mesh src/main.cpp)
target_link_libraries(pqss-proxy-mesh PRIVATE pqss_proxy_mesh)

add_executable(pqss-primitive-mesh-analyze tools/primitive_mesh_analyze.cpp)
target_link_libraries(pqss-primitive-mesh-analyze PRIVATE pqss_proxy_mesh)

if(PQSS_PROXY_MESH_BUILD_REAL_SCENE_COMPARE)
    if(NOT EXISTS "${PQSS_SOURCE_DIR}/src/pqss/include/pqss/model_pool.hpp")
        message(FATAL_ERROR "PQSS source tree not found at PQSS_SOURCE_DIR=${PQSS_SOURCE_DIR}")
    endif()

    add_library(pqss_cpu STATIC
        src/independent_primitive_pool.cpp
        src/primitive_bvh_pool.cpp
        src/pqss_proxy_build.cpp
        src/pqss_proxy_model_pool.cpp
        src/proxy_bvh_control.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/build_cache_path.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/model.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/query_result.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/rect_dist.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/timer.cpp
        ${PQSS_SOURCE_DIR}/src/pqss/src/tri_dist.cpp
    )
    target_include_directories(pqss_cpu PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${PQSS_SOURCE_DIR}/src/pqss/include
    )
    target_compile_features(pqss_cpu PUBLIC cxx_std_20)
    target_compile_definitions(pqss_cpu PRIVATE
        PQSS_CACHE_DIR="${CMAKE_CURRENT_BINARY_DIR}/pqss_cache"
    )

    add_executable(pqss-real-scene-compare tools/real_scene_compare.cpp)
    target_link_libraries(pqss-real-scene-compare PRIVATE pqss_cpu)
    target_compile_features(pqss-real-scene-compare PRIVATE cxx_std_20)

    add_executable(pqss-independent-real-scene-compare tools/independent_real_scene_compare.cpp)
    target_link_libraries(pqss-independent-real-scene-compare PRIVATE pqss_cpu)
    target_compile_features(pqss-independent-real-scene-compare PRIVATE cxx_std_20)

    add_executable(pqss-primitive-bvh-real-scene-compare tools/primitive_bvh_real_scene_compare.cpp)
    target_link_libraries(pqss-primitive-bvh-real-scene-compare PRIVATE pqss_cpu)
    target_compile_features(pqss-primitive-bvh-real-scene-compare PRIVATE cxx_std_20)

    if(MSVC)
        target_compile_options(pqss_cpu PRIVATE /W4 /permissive- /utf-8)
        target_compile_options(pqss-real-scene-compare PRIVATE /W4 /permissive- /utf-8)
        target_compile_options(pqss-independent-real-scene-compare PRIVATE /W4 /permissive- /utf-8)
        target_compile_options(pqss-primitive-bvh-real-scene-compare PRIVATE /W4 /permissive- /utf-8)
    else()
        target_compile_options(pqss_cpu PRIVATE -Wall -Wextra -Wpedantic)
        target_compile_options(pqss-real-scene-compare PRIVATE -Wall -Wextra -Wpedantic)
        target_compile_options(pqss-independent-real-scene-compare PRIVATE -Wall -Wextra -Wpedantic)
        target_compile_options(pqss-primitive-bvh-real-scene-compare PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endif()

if(PQSS_PROXY_MESH_BUILD_TESTS)
    enable_testing()
    add_executable(pqss_proxy_mesh_options_test tests/options_test.cpp)
    target_link_libraries(pqss_proxy_mesh_options_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_proxy_mesh_options_test COMMAND pqss_proxy_mesh_options_test)

    add_executable(pqss_primitive_mesh_analyzer_test tests/primitive_mesh_analyzer_test.cpp)
    target_link_libraries(pqss_primitive_mesh_analyzer_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_primitive_mesh_analyzer_test COMMAND pqss_primitive_mesh_analyzer_test)

    if(PQSS_PROXY_MESH_BUILD_REAL_SCENE_COMPARE)
        add_executable(pqss_independent_primitive_pool_test tests/independent_primitive_pool_test.cpp)
        target_link_libraries(pqss_independent_primitive_pool_test PRIVATE pqss_cpu)
        target_compile_features(pqss_independent_primitive_pool_test PRIVATE cxx_std_20)
        add_test(NAME pqss_independent_primitive_pool_test COMMAND pqss_independent_primitive_pool_test)

        add_executable(pqss_primitive_bvh_pool_test tests/primitive_bvh_pool_test.cpp)
        target_link_libraries(pqss_primitive_bvh_pool_test PRIVATE pqss_cpu)
        target_compile_features(pqss_primitive_bvh_pool_test PRIVATE cxx_std_20)
        add_test(NAME pqss_primitive_bvh_pool_test COMMAND pqss_primitive_bvh_pool_test)
    endif()

    find_package(Python3 COMPONENTS Interpreter)
    if(Python3_Interpreter_FOUND)
        add_test(
            NAME pqss_proxy_mesh_python_test
            COMMAND ${Python3_EXECUTABLE} -m unittest discover -s tests/python -v
        )
        set_tests_properties(pqss_proxy_mesh_python_test PROPERTIES
            ENVIRONMENT "PYTHONPATH=${CMAKE_CURRENT_SOURCE_DIR}/python"
            WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
        )
    endif()
endif()
