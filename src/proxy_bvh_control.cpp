#include "pqss_proxy_mesh/proxy_bvh_build.hpp"

#include <atomic>

namespace pqss_proxy_mesh
{
namespace
{
std::atomic_bool g_proxy_bvh_build_enabled = false;
}

void setProxyBvhBuildEnabled(const bool enabled)
{
    g_proxy_bvh_build_enabled.store(enabled, std::memory_order_relaxed);
}

bool proxyBvhBuildEnabled()
{
    return g_proxy_bvh_build_enabled.load(std::memory_order_relaxed);
}

} // namespace pqss_proxy_mesh
