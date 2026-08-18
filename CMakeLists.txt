cmake_minimum_required(VERSION 3.24)

project(PQSSProxyMesh VERSION 0.1.0 LANGUAGES CXX)

option(PQSS_PROXY_MESH_BUILD_TESTS "Build PQSSProxyMesh tests" ON)
option(PQSS_PROXY_MESH_BUILD_REAL_SCENE_COMPARE "Build the PQSS real-scene comparison tool" ON)
set(PQSS_SOURCE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/../PQSS" CACHE PATH "Path to the PQSS source tree")

add_library(pqss_proxy_mesh STATIC
    src/mesh_pool.cpp
    src/options.cpp
    src/phase0_open_hole.cpp
    src/topology_fill.cpp
    src/halfedge_validation.cpp
    src/phase2_halfedge.cpp
    src/hausdorff_simplifier.cpp
    src/halfedge_surface_analysis.cpp
    src/primitive_mesh_analyzer.cpp
    deps/quickhull/QuickHull.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.engine.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.offset.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.rectclip.cpp
    deps/clipper2/CPP/Clipper2Lib/src/clipper.triangulation.cpp
)

target_include_directories(pqss_proxy_mesh
    PUBLIC
        ${CMAKE_CURRENT_SOURCE_DIR}/include
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/clipper2/CPP/Clipper2Lib/include
        ${CMAKE_CURRENT_SOURCE_DIR}/deps/quickhull
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

add_executable(pqss-topology-fill tools/topology_fill_generate.cpp)
target_link_libraries(pqss-topology-fill PRIVATE pqss_proxy_mesh)

add_executable(pqss-canonicalize-soup tools/canonicalize_triangle_soup.cpp)
target_link_libraries(pqss-canonicalize-soup PRIVATE pqss_proxy_mesh)

add_executable(pqss-planar-hole-cap-preview tools/planar_hole_cap_preview.cpp)
target_link_libraries(pqss-planar-hole-cap-preview PRIVATE pqss_proxy_mesh)

add_executable(pqss-phase0-cap-open-holes tools/phase0_cap_open_holes.cpp)
target_link_libraries(pqss-phase0-cap-open-holes PRIVATE pqss_proxy_mesh)

add_executable(pqss-phase1-from-phase0-batch tools/phase1_from_phase0_batch.cpp)
target_link_libraries(pqss-phase1-from-phase0-batch PRIVATE pqss_proxy_mesh)

add_executable(pqss-hausdorff-simplify tools/hausdorff_simplify.cpp)
target_link_libraries(pqss-hausdorff-simplify PRIVATE pqss_proxy_mesh)

add_executable(pqss-halfedge-surface-analyze tools/halfedge_surface_analyze.cpp)
target_link_libraries(pqss-halfedge-surface-analyze PRIVATE pqss_proxy_mesh)

add_executable(pqss-certify-directed-hausdorff tools/certify_directed_hausdorff.cpp)
target_link_libraries(pqss-certify-directed-hausdorff PRIVATE pqss_proxy_mesh)

add_executable(pqss-refine-proxy-planar tools/refine_proxy_planar.cpp)
target_link_libraries(pqss-refine-proxy-planar PRIVATE pqss_proxy_mesh)

add_executable(pqss-obj-to-halfedge tools/obj_to_halfedge.cpp)
target_link_libraries(pqss-obj-to-halfedge PRIVATE pqss_proxy_mesh)

add_executable(pqss-build-halfedge-surface-manifest tools/build_halfedge_surface_manifest.cpp)
target_compile_features(pqss-build-halfedge-surface-manifest PRIVATE cxx_std_20)

add_executable(pqss-validate-halfedge tools/validate_halfedge.cpp)
target_link_libraries(pqss-validate-halfedge PRIVATE pqss_proxy_mesh)

add_executable(pqss-inspect-halfedge-vertex tools/inspect_halfedge_vertex.cpp)
target_link_libraries(pqss-inspect-halfedge-vertex PRIVATE pqss_proxy_mesh)

add_executable(pqss-inspect-halfedge-face tools/inspect_halfedge_face.cpp)
target_link_libraries(pqss-inspect-halfedge-face PRIVATE pqss_proxy_mesh)

add_executable(pqss-phase2-repair-halfedge tools/phase2_repair_halfedge.cpp)
target_link_libraries(pqss-phase2-repair-halfedge PRIVATE pqss_proxy_mesh)

add_executable(pqss-phase2-concavity-diagnose tools/phase2_concavity_diagnose.cpp)
target_link_libraries(pqss-phase2-concavity-diagnose PRIVATE pqss_proxy_mesh)

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

    add_executable(pqss_topology_fill_test tests/topology_fill_test.cpp)
    target_link_libraries(pqss_topology_fill_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_topology_fill_test COMMAND pqss_topology_fill_test)

    add_executable(pqss_phase0_open_hole_test tests/phase0_open_hole_test.cpp)
    target_link_libraries(pqss_phase0_open_hole_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_phase0_open_hole_test COMMAND pqss_phase0_open_hole_test)

    add_executable(pqss_hausdorff_simplifier_test tests/hausdorff_simplifier_test.cpp)
    target_link_libraries(pqss_hausdorff_simplifier_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_hausdorff_simplifier_test COMMAND pqss_hausdorff_simplifier_test)

    add_executable(pqss_halfedge_surface_analysis_test tests/halfedge_surface_analysis_test.cpp)
    target_link_libraries(pqss_halfedge_surface_analysis_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_halfedge_surface_analysis_test COMMAND pqss_halfedge_surface_analysis_test)

    add_executable(pqss_phase2_halfedge_test tests/phase2_halfedge_test.cpp)
    target_link_libraries(pqss_phase2_halfedge_test PRIVATE pqss_proxy_mesh)
    add_test(NAME pqss_phase2_halfedge_test COMMAND pqss_phase2_halfedge_test)

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

endif()
