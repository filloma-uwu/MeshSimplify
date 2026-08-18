#include <algorithm>
#include <cctype>
#include <exception>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{

std::string readText(const std::filesystem::path& path)
{
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("failed to read " + path.string());
    return {std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>()};
}

void writeText(const std::filesystem::path& path, const std::string& text)
{
    std::ofstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("failed to write " + path.string());
    stream << text;
}

std::string jsonPath(const std::filesystem::path& path)
{
    std::string result = path.generic_string();
    std::string escaped;
    escaped.reserve(result.size());
    for (const char character : result)
    {
        if (character == '\\' || character == '"') escaped.push_back('\\');
        escaped.push_back(character);
    }
    return escaped;
}

std::string numericField(const std::string& json, const std::string& name)
{
    const std::string key = "\"" + name + "\":";
    const std::size_t begin = json.find(key);
    if (begin == std::string::npos)
        throw std::runtime_error("missing numeric field: " + name);
    const std::size_t value_begin = begin + key.size();
    const std::size_t value_end = json.find_first_of(",}", value_begin);
    if (value_end == std::string::npos)
        throw std::runtime_error("unterminated numeric field: " + name);
    return json.substr(value_begin, value_end - value_begin);
}

} // namespace

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path phase1_root;
        std::filesystem::path phase2_root;
        double maximum_directed_hausdorff = 0.0;
        bool has_limit = false;
        for (int index = 1; index < argc; ++index)
        {
            const std::string argument = argv[index];
            const auto value = [&]() -> std::string
            {
                if (++index >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[index];
            };
            if (argument == "--phase1-root") phase1_root = value();
            else if (argument == "--phase2-root") phase2_root = value();
            else if (argument == "--maximum-directed-hausdorff")
            {
                maximum_directed_hausdorff = std::stod(value());
                has_limit = true;
            }
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (phase1_root.empty() || phase2_root.empty() || !has_limit)
            throw std::invalid_argument(
                "usage: pqss-build-halfedge-surface-manifest --phase1-root DIR "
                "--phase2-root DIR --maximum-directed-hausdorff VALUE");

        phase1_root = std::filesystem::absolute(phase1_root);
        phase2_root = std::filesystem::absolute(phase2_root);
        std::vector<int> model_ids;
        for (const auto& entry : std::filesystem::directory_iterator(phase2_root / "models"))
        {
            if (!entry.is_directory()) continue;
            const std::string name = entry.path().filename().string();
            if (name.empty() || !std::all_of(name.begin(), name.end(),
                    [](const unsigned char character) { return std::isdigit(character) != 0; }))
                continue;
            model_ids.push_back(std::stoi(name));
        }
        std::sort(model_ids.begin(), model_ids.end());
        if (model_ids.empty()) throw std::runtime_error("phase2 root has no model outputs");

        for (const int model_id : model_ids)
        {
            const std::filesystem::path model_directory =
                phase2_root / "models" / std::to_string(model_id);
            const std::filesystem::path model_path = model_directory / "model.json";
            if (!std::filesystem::is_regular_file(model_directory / "proxy.obj"))
                throw std::runtime_error("model has no proxy.obj: " + std::to_string(model_id));
            std::string model = readText(model_path);
            while (!model.empty() && std::isspace(static_cast<unsigned char>(model.back())))
                model.pop_back();
            if (model.empty() || model.back() != '}')
                throw std::runtime_error("invalid generated model.json: " + model_path.string());
            const std::filesystem::path phase1_model =
                phase1_root / "models" / std::to_string(model_id);
            if (model.find("\"primitive_count\"") == std::string::npos)
            {
                const std::string phase1_metadata = readText(phase1_model / "model.json");
                const std::string additions =
                    ",\"source_triangles\":" + numericField(phase1_metadata, "source_triangles") +
                    ",\"primitive_count\":" + numericField(model, "regions") +
                    ",\"primitive_types\":{\"polygon\":" +
                        numericField(model, "support_plane_polygons") +
                        ",\"triangle\":" + numericField(model, "fallback_triangles") + "}" +
                    ",\"timings_seconds\":{\"total\":" +
                        numericField(model, "elapsed_seconds") + "}" +
                    ",\"simplification_error\":{\"maximum_distance_limit\":" +
                        numericField(model, "maximum_directed_hausdorff") +
                        ",\"maximum_distance\":" +
                        numericField(model, "certified_directed_hausdorff_upper_bound") + "}";
                std::size_t stats_end = model.find("},\n  \"source\"");
                if (stats_end == std::string::npos) stats_end = model.rfind("}}");
                if (stats_end == std::string::npos)
                    throw std::runtime_error("cannot locate stats object: " + model_path.string());
                model.insert(stats_end, additions);
            }
            if (model.find("\"viewer_stages\"") != std::string::npos)
            {
                writeText(model_path, model + "\n");
                continue;
            }
            model.pop_back();
            const auto relative = [&](const std::filesystem::path& path)
            {
                return jsonPath(std::filesystem::relative(path, model_directory));
            };
            model += ",\n  \"source\": \"" + relative(phase1_model / "source.obj") + "\",\n";
            model += "  \"phase1_halfedge\": \"" +
                     relative(phase1_model / "phase1_halfedge.bin") + "\",\n";
            model += "  \"phase2_recognized_surfaces\": \"proxy.obj\",\n";
            model += "  \"phase4_triangulated\": \"proxy.obj\",\n";
            model += "  \"proxy\": \"proxy.obj\",\n";
            model += "  \"proxy_components\": [],\n";
            model += "  \"viewer_stages\": [\"source\", \"phase1\", \"phase2\", "
                     "\"phase4\", \"split\"]\n}\n";
            writeText(model_path, model);
        }

        std::ofstream manifest(phase2_root / "viewer_manifest.json");
        if (!manifest) throw std::runtime_error("failed to create viewer manifest");
        manifest.precision(17);
        manifest << "{\n  \"algorithm\": \"HalfedgeAdjacentSupportPlaneFixedPointV2\",\n"
                 << "  \"complete\": true,\n  \"model_count\": " << model_ids.size()
                 << ",\n  \"models\": [\n";
        for (std::size_t index = 0; index < model_ids.size(); ++index)
        {
            manifest << "    {\"id\": " << model_ids[index] << ", \"metadata\": "
                     << "\"models/" << model_ids[index] << "/model.json\"}"
                     << (index + 1 == model_ids.size() ? "\n" : ",\n");
        }
        manifest << "  ],\n  \"options\": {\n"
                 << "    \"maximum_directed_hausdorff\": "
                 << maximum_directed_hausdorff << ",\n"
                 << "    \"direction\": \"proxy_to_phase1\",\n"
                 << "    \"phase1_input\": \"" << jsonPath(phase1_root.filename()) << "\",\n"
                 << "    \"process_memory_limit_mib\": 2048,\n"
                 << "    \"external_abort_working_set_mib\": 500,\n"
                 << "    \"maximum_runtime_seconds\": 300,\n"
                 << "    \"no_progress_timeout_seconds\": 120\n"
                 << "  }\n}\n";
        std::cout << "models=" << model_ids.size() << " manifest="
                  << (phase2_root / "viewer_manifest.json").string() << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
