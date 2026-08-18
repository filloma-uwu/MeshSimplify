#include "pqss_proxy_mesh/halfedge_surface_analysis.hpp"
#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#include <psapi.h>
#pragma comment(lib, "psapi.lib")
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

double workingSetMiB()
{
    PROCESS_MEMORY_COUNTERS counters{};
    counters.cb = sizeof(counters);
    return GetProcessMemoryInfo(GetCurrentProcess(), &counters, sizeof(counters))
        ? static_cast<double>(counters.WorkingSetSize) / (1024.0 * 1024.0)
        : -1.0;
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
        std::cerr << "monitor: memory_limit_mib=2048 stage=started\n";
        std::filesystem::path input;
        std::filesystem::path output;
        pqss_proxy_mesh::HalfedgeSurfaceAnalysisOptions options;
        bool has_limit = false;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--phase1-halfedge") input = value();
            else if (argument == "--output-dir") output = value();
            else if (argument == "--maximum-directed-hausdorff")
            {
                options.maximum_directed_hausdorff = std::stod(value());
                has_limit = true;
            }
            else
                throw std::invalid_argument("unknown argument: " + argument);
        }
        if (input.empty() || output.empty() || !has_limit)
            throw std::invalid_argument(
                "usage: pqss-halfedge-surface-analyze --phase1-halfedge FILE "
                "--output-dir DIR --maximum-directed-hausdorff VALUE");
        const pqss_proxy_mesh::OrientedSurfaceMesh mesh =
            pqss_proxy_mesh::readAnalysisHalfedgeMesh(input);
#ifdef _WIN32
        std::cerr << "monitor: stage=phase1_loaded working_set_mib="
                  << workingSetMiB() << '\n';
#endif
        const auto stats = pqss_proxy_mesh::analyzeHalfedgeSurface(mesh, output, options);
#ifdef _WIN32
        std::cerr << "monitor: stage=analysis_finished working_set_mib="
                  << workingSetMiB() << '\n';
#endif
        std::cout << "phase1_triangles=" << stats.phase1_triangles << '\n'
                  << "proxy_triangles=" << stats.final_triangles << '\n'
                  << "planar_regions=" << stats.planar_regions << '\n'
                  << "emitted_planar_polygons=" << stats.emitted_planar_polygons << '\n'
                   << "fallback_triangles=" << stats.fallback_triangles << '\n'
                   << "coverage_certified=" << (stats.coverage_certified ? "true" : "false")
                   << '\n'
                   << "global_hausdorff_certified="
                   << (stats.global_hausdorff_certified ? "true" : "false") << '\n'
                   << "certified_directed_hausdorff_upper_bound="
                   << stats.certified_directed_hausdorff_upper_bound << '\n'
                   << "elapsed_seconds=" << stats.elapsed_seconds << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
