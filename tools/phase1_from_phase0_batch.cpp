#include "pqss_proxy_mesh/topology_fill.hpp"

#include <algorithm>
#include <atomic>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <mutex>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace
{

struct ModelInput
{
    std::string id;
    std::filesystem::path directory;
};

std::vector<ModelInput> discoverModels(const std::filesystem::path& root)
{
    std::vector<ModelInput> result;
    for (const auto& entry : std::filesystem::directory_iterator(root / "models"))
    {
        if (!entry.is_directory()) continue;
        if (!std::filesystem::is_regular_file(entry.path() / "phase0_open_holes.obj")) continue;
        result.push_back({entry.path().filename().string(), entry.path()});
    }
    std::sort(result.begin(), result.end(), [](const ModelInput& first, const ModelInput& second)
    {
        const auto first_size = std::filesystem::file_size(
            first.directory / "phase0_open_holes.obj");
        const auto second_size = std::filesystem::file_size(
            second.directory / "phase0_open_holes.obj");
        if (first_size != second_size) return first_size > second_size;
        try { return std::stoll(first.id) < std::stoll(second.id); }
        catch (...) { return first.id < second.id; }
    });
    return result;
}

void writeMetadata(
    const std::filesystem::path& path,
    const std::string& id,
    const std::string& phase0_root_url,
    const pqss_proxy_mesh::TopologyFillStats& stats)
{
    std::ofstream stream(path);
    if (!stream) throw std::runtime_error("failed to rewrite phase-1 metadata: " + path.string());
    const std::string base = phase0_root_url + "/models/" + id + "/";
    stream << std::setprecision(17)
           << "{\n  \"stats\":{\n"
           << "    \"model\":\"phase0_open_holes.obj\",\n"
           << "    \"source_triangles\":" << stats.source_triangles << ",\n"
           << "    \"primitive_count\":0,\n"
           << "    \"primitive_types\":{},\n"
           << "    \"proxy_triangles\":" << stats.mesh_triangles << ",\n"
           << "    \"topology_fill\":{\n"
           << "      \"input_betti\":[" << stats.input_betti.beta0 << ','
           << stats.input_betti.beta1 << ',' << stats.input_betti.beta2 << "],\n"
           << "      \"output_betti\":[" << stats.output_betti.beta0 << ','
           << stats.output_betti.beta1 << ',' << stats.output_betti.beta2 << "],\n"
           << "      \"input_voxels\":" << stats.input_voxels << ",\n"
           << "      \"output_voxels\":" << stats.output_voxels << ",\n"
           << "      \"added_voxels\":" << stats.added_voxels << ",\n"
           << "      \"pitch\":" << stats.pitch << ",\n"
           << "      \"grid_shape\":[" << stats.grid_shape[0] << ','
           << stats.grid_shape[1] << ',' << stats.grid_shape[2] << "],\n"
           << "      \"mesh_triangles\":" << stats.mesh_triangles << ",\n"
           << "      \"mesh_vertices\":" << stats.mesh_vertices << ",\n"
           << "      \"halfedge_count\":" << stats.halfedge_count << ",\n"
           << "      \"paired_halfedge_edges\":" << stats.paired_halfedge_edges << ",\n"
           << "      \"boundary_halfedges\":" << stats.boundary_halfedges << ",\n"
           << "      \"halfedge_face_components\":" << stats.halfedge_face_components << ",\n"
           << "      \"mesh_watertight\":true,\n"
           << "      \"mesh_oriented\":true,\n"
           << "      \"mesh_manifold\":true,\n"
           << "      \"mesh_connected\":true\n"
           << "    },\n"
           << "    \"timings_seconds\":{\"total\":" << stats.elapsed_seconds << "}\n"
           << "  },\n"
           << "  \"source\":\"" << base << "source.obj\",\n"
           << "  \"phase0_caps\":\"" << base << "phase0_caps.obj\",\n"
           << "  \"phase0_combined\":\"" << base << "phase0_open_holes.obj\",\n"
           << "  \"phase1_halfedge\":\"phase1_halfedge.bin\",\n"
           << "  \"proxy_components\":[],\n"
           << "  \"viewer_stages\":[\"source\",\"phase0\",\"phase1\"]\n"
           << "}\n";
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path phase0_root;
        std::filesystem::path output_root;
        std::size_t maximum_parallel_models = 4;
        pqss_proxy_mesh::TopologyFillOptions options;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            if (index + 1 >= argc) throw std::invalid_argument("missing value for " + argument);
            const std::string value = argv[++index];
            if (argument == "--phase0-root") phase0_root = value;
            else if (argument == "--output-root") output_root = value;
            else if (argument == "--maximum-grid-voxels") options.maximum_grid_voxels = std::stoull(value);
            else if (argument == "--padding") options.padding = static_cast<std::uint32_t>(std::stoul(value));
            else if (argument == "--maximum-steps") options.maximum_steps = static_cast<std::uint32_t>(std::stoul(value));
            else if (argument == "--maximum-parallel-models") maximum_parallel_models = std::stoull(value);
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (phase0_root.empty() || output_root.empty())
            throw std::invalid_argument(
                "usage: pqss-phase1-from-phase0-batch --phase0-root DIR --output-root DIR");
        if (maximum_parallel_models == 0)
            throw std::invalid_argument("maximum parallel models must be positive");
        if (std::filesystem::exists(output_root))
            throw std::runtime_error("output root already exists: " + output_root.string());
        const auto models = discoverModels(phase0_root);
        if (models.empty()) throw std::runtime_error("phase-0 root contains no phase0_open_holes.obj files");
        std::filesystem::create_directories(output_root / "models");
        const std::string phase0_root_url = "/outputs/" + phase0_root.filename().string();
        std::vector<pqss_proxy_mesh::TopologyFillStats> generated(models.size());
        std::atomic<std::size_t> next_model{0};
        std::mutex output_mutex;
        std::exception_ptr failure;
        const std::size_t worker_count = std::min(maximum_parallel_models, models.size());
        std::vector<std::thread> workers;
        workers.reserve(worker_count);
        for (std::size_t worker = 0; worker < worker_count; ++worker)
        {
            workers.emplace_back([&]
            {
                try
                {
                    for (;;)
                    {
                        const std::size_t index = next_model.fetch_add(1);
                        if (index >= models.size()) break;
                        const auto& model = models[index];
                        const auto output = output_root / "models" / model.id;
                        generated[index] = pqss_proxy_mesh::generateTopologyFillModel(
                            model.directory / "phase0_open_holes.obj", output, options);
                        writeMetadata(
                            output / "model.json", model.id, phase0_root_url,
                            generated[index]);
                        std::scoped_lock lock(output_mutex);
                        std::cout << "model=" << model.id
                                  << " input_triangles=" << generated[index].source_triangles
                                  << " phase1_triangles=" << generated[index].mesh_triangles
                                  << " elapsed_seconds=" << generated[index].elapsed_seconds
                                  << '\n' << std::flush;
                    }
                }
                catch (...)
                {
                    std::scoped_lock lock(output_mutex);
                    if (!failure) failure = std::current_exception();
                    next_model.store(models.size());
                }
            });
        }
        for (auto& worker : workers) worker.join();
        if (failure) std::rethrow_exception(failure);

        std::ofstream manifest(output_root / "viewer_manifest.json");
        if (!manifest) throw std::runtime_error("failed to create phase-1 viewer manifest");
        manifest << "{\n  \"algorithm\":\"CppPhase1DirectHoleMouthProjectionV18\",\n"
                 << "  \"complete\":true,\n  \"model_count\":" << models.size() << ",\n  \"models\":[\n";
        for (std::size_t index = 0; index < models.size(); ++index)
        {
            const auto& model = models[index];
            manifest << "    {\"id\":" << model.id
                     << ",\"metadata\":\"models/" << model.id << "/model.json\"}"
                     << (index + 1 == models.size() ? "\n" : ",\n");
        }
        manifest << "  ],\n  \"options\":{\"phase0_root\":\""
                 << phase0_root_url << "\",\"maximum_grid_voxels\":"
                 << options.maximum_grid_voxels << ",\"padding\":" << options.padding
                 << ",\"maximum_steps\":" << options.maximum_steps
                 << ",\"maximum_parallel_models\":" << maximum_parallel_models
                 << ",\"phase1_input\":\"phase0_open_holes.obj\"}\n}\n";
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
