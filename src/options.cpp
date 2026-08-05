#include "pqss_proxy_mesh/options.hpp"

#include <cmath>

namespace pqss_proxy_mesh
{

std::optional<std::string> validateOptions(const SimplificationOptions& options)
{
    if (options.max_pqss_bvh_depth == 0)
    {
        return "max_pqss_bvh_depth must be greater than zero";
    }
    if (!std::isfinite(options.relative_alpha) || options.relative_alpha <= 0.0)
    {
        return "relative_alpha must be finite and greater than zero";
    }
    if (!std::isfinite(options.relative_offset) || options.relative_offset <= 0.0)
    {
        return "relative_offset must be finite and greater than zero";
    }
    return std::nullopt;
}

} // namespace pqss_proxy_mesh
