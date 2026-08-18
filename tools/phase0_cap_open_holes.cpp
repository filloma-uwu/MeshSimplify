#include "pqss_proxy_mesh/phase0_open_hole.hpp"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <algorithm>
#include <stdexcept>
#include <string>
#include <vector>

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path input;
        std::filesystem::path source_root;
        std::filesystem::path output;
        std::filesystem::path output_root;
        std::string model_id;
        pqss_proxy_mesh::Phase0OpenHoleOptions options;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc)
                    throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--input") input = value();
            else if (argument == "--source-root") source_root = value();
            else if (argument == "--output-dir") output = value();
            else if (argument == "--output-root") output_root = value();
            else if (argument == "--model-id") model_id = value();
            else if (argument == "--projection-relative-tolerance")
                options.projection_relative_tolerance = std::stod(value());
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (!output_root.empty())
        {
            if (model_id.empty() && source_root.empty())
                throw std::invalid_argument("--output-root requires --model-id");
            if (!model_id.empty()) output = output_root / "models" / model_id;
        }
        if (!source_root.empty())
        {
            if (output_root.empty())
                throw std::invalid_argument("--source-root requires --output-root");
            struct BatchModel { std::string id; std::filesystem::path source; };
            std::vector<BatchModel> models;
            for (const auto& entry : std::filesystem::directory_iterator(source_root))
            {
                if (!entry.is_directory()) continue;
                const auto candidate = entry.path() / "source.obj";
                if (!std::filesystem::is_regular_file(candidate)) continue;
                models.push_back({entry.path().filename().string(), candidate});
            }
            std::sort(models.begin(), models.end(), [](const BatchModel& first,
                                                       const BatchModel& second)
            {
                try { return std::stoll(first.id) < std::stoll(second.id); }
                catch (...) { return first.id < second.id; }
            });
            if (models.empty()) throw std::invalid_argument("source root contains no model source.obj files");
            std::filesystem::create_directories(output_root);
            std::ofstream manifest(output_root / "viewer_manifest.json");
            if (!manifest) throw std::runtime_error("failed to create phase-0 viewer manifest");
            manifest << "{\n  \"algorithm\":\"CppOriginalObjOpenBoundaryCapsV2\",\n"
                     << "  \"complete\":false,\n  \"model_count\":" << models.size() << ",\n  \"models\":[\n";
            for (std::size_t index = 0; index < models.size(); ++index)
            {
                const auto& model = models[index];
                const auto model_output = output_root / "models" / model.id;
                const auto stats = pqss_proxy_mesh::generatePhase0OpenHoleModel(
                    model.source, model_output, options);
                std::cout << "model=" << model.id
                          << " source_triangles=" << stats.source_triangles
                          << " reality_boundary_edges=" << stats.reality_boundary_edges
                          << " closed_boundary_loops=" << stats.closed_boundary_loops
                          << " capped_loops=" << stats.capped_loops
                          << " cap_triangles=" << stats.cap_triangles
                          << " elapsed_seconds=" << stats.elapsed_seconds << '\n';
                manifest << "    {\"id\":";
                try { manifest << std::stoll(model.id); }
                catch (...) { manifest << '"' << model.id << '"'; }
                manifest << ",\"metadata\":\"models/" << model.id
                         << "/model.json\"}"
                         << (index + 1 == models.size() ? "\n" : ",\n");
            }
            manifest << "  ],\n  \"options\":{\"input\":\"original_obj\","
                     << "\"occupancy_used\":false,\"projection_relative_tolerance\":"
                     << std::setprecision(17) << options.projection_relative_tolerance
                     << "},\n  \"complete\":true\n}\n";
            return 0;
        }
        if (input.empty() || output.empty())
            throw std::invalid_argument(
                "usage: pqss-phase0-cap-open-holes --input source.obj "
                "--output-dir DIR [--projection-relative-tolerance VALUE]");
        const auto stats = pqss_proxy_mesh::generatePhase0OpenHoleModel(
            input, output, options);
        if (!output_root.empty())
        {
            std::filesystem::create_directories(output_root);
            std::ofstream manifest(output_root / "viewer_manifest.json");
            if (!manifest) throw std::runtime_error("failed to create phase-0 viewer manifest");
            manifest << std::setprecision(17)
                     << "{\n  \"algorithm\":\"CppOriginalObjOpenBoundaryCapsV2\",\n"
                     << "  \"complete\":true,\n  \"model_count\":1,\n"
                     << "  \"models\":[{\"id\":" << model_id
                     << ",\"metadata\":\"models/" << model_id
                     << "/model.json\"}],\n"
                     << "  \"options\":{\"input\":\"original_obj\","
                     << "\"occupancy_used\":false,"
                     << "\"projection_relative_tolerance\":"
                     << options.projection_relative_tolerance << "}\n}\n";
        }
        std::cout << "source_triangles=" << stats.source_triangles << '\n'
                  << "reality_boundary_edges=" << stats.reality_boundary_edges << '\n'
                  << "closed_boundary_loops=" << stats.closed_boundary_loops << '\n'
                  << "capped_loops=" << stats.capped_loops << '\n'
                  << "cap_triangles=" << stats.cap_triangles << '\n'
                  << "elapsed_seconds=" << stats.elapsed_seconds << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
