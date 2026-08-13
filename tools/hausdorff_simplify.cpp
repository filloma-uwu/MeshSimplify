#include "pqss_proxy_mesh/hausdorff_simplifier.hpp"

#include <cmath>
#include <exception>
#include <filesystem>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv)
{
    try
    {
        std::filesystem::path input;
        std::filesystem::path source;
        std::filesystem::path output;
        pqss_proxy_mesh::HausdorffSimplificationOptions options;
        bool has_limit = false;
        for (int i = 1; i < argc; ++i)
        {
            const std::string argument = argv[i];
            const auto value = [&]() -> std::string
            {
                if (++i >= argc) throw std::invalid_argument("missing value for " + argument);
                return argv[i];
            };
            if (argument == "--phase1-halfedge") input = value();
            else if (argument == "--source") source = value();
            else if (argument == "--output-dir") output = value();
            else if (argument == "--model-id") options.model_id = value();
            else if (argument == "--maximum-directed-hausdorff")
            {
                options.maximum_directed_hausdorff = std::stod(value());
                has_limit = true;
            }
            else if (argument == "--maximum-certificate-depth")
                options.maximum_certificate_depth = static_cast<std::uint32_t>(std::stoul(value()));
            else if (argument == "--no-coplanar-union") options.enable_exact_coplanar_union = false;
            else if (argument == "--no-convex-hull") options.enable_convex_hull = false;
            else if (argument == "--no-adaptive-convex-cover")
                options.enable_adaptive_convex_cover = false;
            else if (argument == "--no-kdop") options.enable_discrete_orientation_polytopes = false;
            else if (argument == "--no-box") options.enable_axis_aligned_box = false;
            else throw std::invalid_argument("unknown argument: " + argument);
        }
        if (input.empty() || source.empty() || output.empty() || !has_limit)
            throw std::invalid_argument(
                "usage: pqss-hausdorff-simplify --phase1-halfedge FILE --source FILE "
                "--output-dir DIR --maximum-directed-hausdorff VALUE");
        const auto stats = pqss_proxy_mesh::simplifyPhase1Halfedge(input, source, output, options);
        std::cout << "phase1_triangles=" << stats.phase1_triangles << '\n'
                  << "proxy_triangles=" << stats.final_triangles << '\n'
                  << "selected_candidate=" << stats.selected_candidate << '\n'
                  << "elapsed_seconds=" << stats.elapsed_seconds << '\n';
        return 0;
    }
    catch (const std::exception& error)
    {
        std::cerr << error.what() << '\n';
        return 1;
    }
}
