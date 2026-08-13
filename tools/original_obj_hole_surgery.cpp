#include "pqss_proxy_mesh/topology_fill.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace
{
#ifdef _WIN32
HANDLE installMemoryLimit()
{
    constexpr SIZE_T memory_limit = SIZE_T{2} * 1024 * 1024 * 1024;
    HANDLE job = CreateJobObjectW(nullptr, nullptr);
    if (!job) throw std::runtime_error("failed to create memory-limit job");
    JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
    limits.BasicLimitInformation.LimitFlags = JOB_OBJECT_LIMIT_PROCESS_MEMORY;
    limits.ProcessMemoryLimit = memory_limit;
    if (!SetInformationJobObject(job, JobObjectExtendedLimitInformation,
            &limits, sizeof(limits)) ||
        !AssignProcessToJobObject(job, GetCurrentProcess()))
    {
        CloseHandle(job);
        throw std::runtime_error("failed to install 2 GiB process memory limit");
    }
    return job;
}
#endif
}

int main(int argc, char** argv)
{
    try
    {
#ifdef _WIN32
        const HANDLE memory_job = installMemoryLimit();
#endif
        if (argc != 5)
            throw std::invalid_argument(
                "usage: pqss-original-obj-hole-surgery source.obj filled.vox cavity-labels.vox output.halfedge");
        const auto source = pqss_proxy_mesh::readTriangleSoupObj(argv[1]);
        const auto filled = pqss_proxy_mesh::readVoxelGrid(argv[2]);
        const auto cavity_labels = pqss_proxy_mesh::readVoxelGrid(argv[3]);
        pqss_proxy_mesh::TriangleVoxelProvenance provenance;
        const auto original = pqss_proxy_mesh::voxelizeTriangleSoupOnGrid(
            source, filled, &provenance);
        std::cerr << "monitor: memory_limit_mib=2048 stage=loaded\n";
        pqss_proxy_mesh::OriginalObjHoleSurgeryStats stats;
        const auto surgery_soup = pqss_proxy_mesh::buildOriginalObjHoleSurgery(
            filled, cavity_labels, original, provenance, source, &stats);
        const auto surgery = pqss_proxy_mesh::canonicalizeCoplanarTriangleSoup(
            surgery_soup);
        pqss_proxy_mesh::AnalysisHalfedgeStats halfedge_stats;
        const auto result = pqss_proxy_mesh::buildAnalysisHalfedgeMesh(
            surgery, &halfedge_stats);
        pqss_proxy_mesh::writeAnalysisHalfedgeMesh(argv[4], result);
        std::cout << "{\"source_voxels\":" << stats.source_voxels
                  << ",\"filled_voxels\":" << stats.filled_voxels
                  << ",\"added_voxels\":" << stats.added_voxels
                  << ",\"source_triangles\":" << stats.source_triangles
                  << ",\"retained_triangles\":" << stats.retained_triangles
                  << ",\"removed_internal_triangles\":"
                  << stats.removed_internal_triangles
                  << ",\"removed_closed_cavity_triangles\":"
                  << stats.removed_closed_cavity_triangles
                  << ",\"removed_closed_cavity_components\":"
                  << stats.removed_closed_cavity_components
                  << ",\"closed_cavity_voxels\":"
                  << stats.closed_cavity_voxels
                  << ",\"ambiguous_triangles\":" << stats.ambiguous_triangles
                  << ",\"interface_edges\":" << stats.interface_edges
                  << ",\"closed_loops\":" << stats.closed_loops
                  << ",\"cap_triangles\":" << stats.cap_triangles << "}\n";
        std::cerr << "monitor: halfedge_components=" << halfedge_stats.face_components
                  << " paired_edges=" << halfedge_stats.paired_edges
                  << " boundary_halfedges=" << halfedge_stats.boundary_halfedges
                  << " duplicate_triangles="
                  << halfedge_stats.dropped_duplicate_triangles << '\n';
#ifdef _WIN32
        CloseHandle(memory_job);
#endif
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
