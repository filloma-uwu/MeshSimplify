#include "pqss_proxy_mesh/phase2_halfedge.hpp"

#include <exception>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path input;
        std::filesystem::path output;
        std::filesystem::path output_root;
        pqss_proxy_mesh::Phase2HalfedgeOptions options;
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
            else if (argument == "--output-root") output_root = value();
            else if (argument == "--model-id") options.model_id = value();
            else if (argument == "--maximum-directed-hausdorff")
                options.maximum_directed_hausdorff = std::stod(value());
            else if (argument == "--maximum-certificate-depth")
                options.maximum_certificate_depth =
                    static_cast<std::uint32_t>(std::stoul(value()));
            else if (argument == "--numerical-tolerance")
                options.numerical_tolerance = std::stod(value());
            else if (argument == "--coplanar-relative-tolerance")
                options.coplanar_relative_tolerance = std::stod(value());
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (!output_root.empty())
        {
            if (options.model_id.empty())
                throw std::invalid_argument("--output-root requires --model-id");
            output = output_root / "models" / options.model_id;
        }
        if (input.empty() || output.empty())
            throw std::invalid_argument(
                "usage: pqss-phase2-repair-halfedge --phase1-halfedge FILE "
                "--output-dir DIR [--model-id ID] "
                "[--maximum-directed-hausdorff VALUE]");
        const auto stats = pqss_proxy_mesh::generatePhase2Halfedge(input, output, options);
        if (!output_root.empty())
        {
            std::filesystem::create_directories(output_root);
            std::ofstream manifest(output_root / "viewer_manifest.json");
            if (!manifest) throw std::runtime_error("failed to create viewer manifest");
            manifest << "{\n  \"algorithm\":\"CppPhase2ConservativeOutwardEnvelopeV1\",\n"
                     << "  \"complete\":true,\n  \"model_count\":1,\n"
                     << "  \"models\":[{\"id\":" << options.model_id
                     << ",\"metadata\":\"models/" << options.model_id
                     << "/model.json\"}],\n"
                     << "  \"options\":{\"phase2\":\"conservative_outward_envelope\","
                     << "\"maximum_directed_hausdorff\":" << std::setprecision(17)
                     << options.maximum_directed_hausdorff
                     << ",\"direction\":\"proxy_to_phase1\"}\n}\n";
        }
        std::cout << "input_triangles=" << stats.input_triangles << '\n'
                  << "output_triangles=" << stats.output_triangles << '\n'
                  << "boundary_halfedges=" << stats.boundary_halfedges << '\n'
                  << "face_components=" << stats.face_components << '\n'
                  << "candidate_faces_considered="
                  << stats.candidate_faces_considered << '\n'
                  << "candidate_faces_accepted="
                  << stats.candidate_faces_accepted << '\n'
                  << "candidate_faces_rejected_disconnected="
                  << stats.candidate_faces_rejected_disconnected << '\n'
                  << "source_faces_retriangulated="
                  << stats.source_faces_retriangulated << '\n'
                  << "selected_candidate=" << stats.selected_candidate << '\n'
                  << "hausdorff_lower_bound="
                  << stats.directed_hausdorff.lower_bound << '\n'
                  << "hausdorff_upper_bound="
                  << stats.directed_hausdorff.upper_bound << '\n'
                  << "elapsed_seconds=" << stats.elapsed_seconds << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
