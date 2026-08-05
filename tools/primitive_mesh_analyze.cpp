#include "pqss_proxy_mesh/primitive_mesh_analyzer.hpp"

#include <cstdlib>
#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <string_view>

#ifdef _WIN32
#define NOMINMAX
#include <windows.h>
#endif

namespace
{

void usage()
{
    std::cerr
        << "usage: pqss-primitive-mesh-analyze --input model.obj --output-dir directory\n"
        << "       [--primitive-types polygon,surface] [--frustum-segments 24]\n"
        << "       [--legacy-analysis-strength 0.5] [--maximum-added-volume-ratio 0.01]\n"
        << "       [--maximum-group-box-added-volume-ratio 0.05]\n"
        << "       [--protrusion-max-area-excess-ratio 0.03] [--diagnostic-volume-envelope]\n"
        << "       [--maximum-open-error-distance value]\n"
        << "       [--maximum-process-memory-gb 2]\n"
        ;
}

void applyProcessMemoryLimit(const double gigabytes)
{
    constexpr double hard_maximum_gigabytes = 2.0;
    if (!std::isfinite(gigabytes) || gigabytes <= 0.0 ||
        gigabytes > hard_maximum_gigabytes)
        throw std::invalid_argument(
            "maximum process memory must be finite, positive, and no greater than 2 GB");
#ifdef _WIN32
    static HANDLE memory_job = nullptr;
    memory_job = CreateJobObjectW(nullptr, nullptr);
    if (!memory_job)
        throw std::runtime_error("failed to create the process-memory job object");

    JOBOBJECT_EXTENDED_LIMIT_INFORMATION limits{};
    limits.BasicLimitInformation.LimitFlags =
        JOB_OBJECT_LIMIT_PROCESS_MEMORY | JOB_OBJECT_LIMIT_KILL_ON_JOB_CLOSE;
    constexpr double bytes_per_gigabyte = 1024.0 * 1024.0 * 1024.0;
    limits.ProcessMemoryLimit = static_cast<SIZE_T>(gigabytes * bytes_per_gigabyte);
    if (!SetInformationJobObject(
            memory_job, JobObjectExtendedLimitInformation,
            &limits, sizeof(limits)) ||
        !AssignProcessToJobObject(memory_job, GetCurrentProcess()))
    {
        const DWORD error = GetLastError();
        CloseHandle(memory_job);
        memory_job = nullptr;
        throw std::runtime_error(
            "failed to install the process-memory limit (Windows error " +
            std::to_string(error) + ")");
    }
#else
    (void)gigabytes;
#endif
}

bool containsType(const std::string& value, const std::string_view type)
{
    std::size_t begin = 0;
    while (begin <= value.size())
    {
        const std::size_t end = value.find(',', begin);
        const std::string_view field(value.data() + begin,
                                     (end == std::string::npos ? value.size() : end) - begin);
        if (field == type) return true;
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    return false;
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path input;
        std::filesystem::path output;
        pqss_proxy_mesh::PrimitiveMeshAnalysisOptions options;
        double maximum_process_memory_gb = 2.0;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--input") input = value();
            else if (argument == "--output-dir") output = value();
            else if (argument == "--primitive-types")
            {
                const std::string types = value();
                options.allow_polygon = containsType(types, "polygon");
                options.allow_frustum = containsType(types, "surface") ||
                                         containsType(types, "frustum");
            }
            else if (argument == "--frustum-segments")
            {
                options.frustum_segments = static_cast<std::uint32_t>(std::stoul(value()));
            }
            else if (argument == "--analysis-strength")
            {
                options.analysis_strength = std::stod(value());
                options.uniform_structure_policy = false;
            }
            else if (argument == "--legacy-analysis-strength")
            {
                options.analysis_strength = std::stod(value());
                options.uniform_structure_policy = false;
            }
            else if (argument == "--maximum-added-volume-ratio")
            {
                options.maximum_added_volume_ratio = std::stod(value());
            }
            else if (argument == "--maximum-group-box-added-volume-ratio")
            {
                options.maximum_group_box_added_volume_ratio = std::stod(value());
            }
            else if (argument == "--protrusion-max-area-excess-ratio")
            {
                options.protrusion_max_area_excess_ratio = std::stod(value());
            }
            else if (argument == "--diagnostic-volume-envelope")
            {
                options.enable_volume_evaluated_envelope = true;
            }
            else if (argument == "--maximum-open-error-distance")
            {
                options.maximum_open_error_distance = std::stod(value());
                if (!std::isfinite(options.maximum_open_error_distance) ||
                    options.maximum_open_error_distance < 0.0)
                    throw std::invalid_argument(
                        "maximum open error distance must be finite and non-negative");
            }
            else if (argument == "--maximum-process-memory-gb")
            {
                maximum_process_memory_gb = std::stod(value());
            }
            else
            {
                usage();
                throw std::invalid_argument("unknown argument: " + argument);
            }
        }
        if (input.empty() || output.empty())
        {
            usage();
            return 2;
        }
        applyProcessMemoryLimit(maximum_process_memory_gb);
        const auto stats = pqss_proxy_mesh::analyzePrimitiveMeshObj(input, output, options);
        std::cout << "source_triangles=" << stats.source_triangles << '\n'
                  << "primitives=" << stats.primitive_count << '\n'
                  << "proxy_triangles=" << stats.proxy_triangles << '\n'
                  << "analysis_seconds=" << stats.analysis_seconds << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
