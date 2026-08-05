#include "pqss_proxy_mesh/options.hpp"

#include <iostream>

int main()
{
    const pqss_proxy_mesh::SimplificationOptions options;
    if (const auto error = pqss_proxy_mesh::validateOptions(options))
    {
        std::cerr << *error << '\n';
        return 1;
    }

    std::cout << "PQSSProxyMesh project initialized\n"
              << "mode=ConservativeOuter\n"
              << "input=model_pool\n"
              << "pqss_subdivision=model_pool_default\n"
              << "max_pqss_bvh_depth=" << options.max_pqss_bvh_depth << '\n';
    return 0;
}
