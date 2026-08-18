#include "pqss_proxy_mesh/topology_fill.hpp"

#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

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

std::vector<int> parseModelIds(const std::string& value)
{
    std::vector<int> result;
    std::size_t begin = 0;
    while (begin < value.size())
    {
        const std::size_t end = value.find(',', begin);
        result.push_back(std::stoi(value.substr(begin, end - begin)));
        if (end == std::string::npos) break;
        begin = end + 1;
    }
    if (result.empty()) throw std::invalid_argument("model ID list is empty");
    return result;
}

void usage()
{
    std::cerr
        << "usage: pqss-topology-fill --source-dir obj_directory --model-ids 2,3,4\n"
        << "       --output-root directory [--maximum-grid-voxels 12000000] [--padding 4]\n"
        << "   or: pqss-topology-fill --input model.obj --output-dir directory\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
#ifdef _WIN32
        const HANDLE memory_job = installMemoryLimit();
#endif
        std::cerr << "monitor: memory_limit_mib=2048 stage=started\n";
        std::filesystem::path source_directory;
        std::filesystem::path output_root;
        std::filesystem::path input;
        std::filesystem::path output_directory;
        std::vector<int> model_ids;
        pqss_proxy_mesh::TopologyFillOptions options;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc)
                    throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--source-dir") source_directory = value();
            else if (argument == "--model-ids") model_ids = parseModelIds(value());
            else if (argument == "--output-root") output_root = value();
            else if (argument == "--input") input = value();
            else if (argument == "--output-dir") output_directory = value();
            else if (argument == "--maximum-grid-voxels")
                options.maximum_grid_voxels = std::stoull(value());
            else if (argument == "--padding")
                options.padding = static_cast<std::uint32_t>(std::stoul(value()));
            else if (argument == "--maximum-steps")
                options.maximum_steps = static_cast<std::uint32_t>(std::stoul(value()));
            else if (argument == "--trace-face")
                options.trace_face = static_cast<std::uint32_t>(std::stoul(value()));
            else
            {
                usage();
                throw std::invalid_argument("unknown argument: " + argument);
            }
        }

        if (!source_directory.empty() || !output_root.empty() || !model_ids.empty())
        {
            if (source_directory.empty() || output_root.empty() || model_ids.empty() ||
                !input.empty() || !output_directory.empty())
            {
                usage();
                return 2;
            }
            const auto results = pqss_proxy_mesh::generateTopologyFillBatch(
                source_directory, model_ids, output_root, options);
            std::size_t triangles = 0;
            for (const auto& result : results) triangles += result.mesh_triangles;
            std::cout << "models=" << results.size() << '\n'
                      << "phase1_mesh_triangles=" << triangles << '\n';
            return 0;
        }

        if (input.empty() || output_directory.empty())
        {
            usage();
            return 2;
        }
        const auto result = pqss_proxy_mesh::generateTopologyFillModel(
            input, output_directory, options);
        std::cout << "input_betti=" << result.input_betti.beta0 << ','
                  << result.input_betti.beta1 << ',' << result.input_betti.beta2 << '\n'
                  << "output_betti=" << result.output_betti.beta0 << ','
                  << result.output_betti.beta1 << ',' << result.output_betti.beta2 << '\n'
                  << "phase1_mesh_triangles=" << result.mesh_triangles << '\n';
        return 0;
    }
    catch (const std::exception& exception)
    {
        std::cerr << exception.what() << '\n';
        return 1;
    }
}
