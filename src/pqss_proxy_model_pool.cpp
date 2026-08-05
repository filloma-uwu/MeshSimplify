#include "pqss/model_pool.hpp"

#include <algorithm>
#include <array>
#include <cctype>
#include <charconv>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <format>
#include <functional>
#include <fstream>
#include <iostream>
#include <limits>
#include <string>
#include <string_view>
#include <system_error>
#include <type_traits>

#include "pqss/build_cache_path.hpp"
#include "pqss/dist.hpp"
#include "pqss/utils/mat_vec.hpp"
#include "pqss/utils/timer.hpp"

namespace pqss
{

namespace
{
using build::BV;
using build::SubTri;
using build::Tri;
} // namespace

// ============================================================================
// Build internal helpers
// ============================================================================

Real ModelPool::TriArea(const Vec3& p1, const Vec3& p2, const Vec3& p3) const
{
    return Vlength_sq(VcrossV(p2 - p1, p3 - p1));
}

Real ModelPool::ComputePoolMeanAreaThreshold() const
{
    Real        all_area     = k_zero;
    std::size_t num_all_tris = 0;
    for (const build::Model& model : m_models)
    {
        for (const Tri& tri : model.OriginalTris())
        {
            all_area += TriArea(tri.p1, tri.p2, tri.p3);
        }
        num_all_tris += model.OriginalTris().size();
    }
    if (num_all_tris == 0)
    {
        return k_zero;
    }
    return all_area / static_cast<Real>(num_all_tris);
}

Real ModelPool::ComputePoolAreaThreshold() const
{
    if (m_area_threshold_override > k_zero)
    {
        return m_area_threshold_override;
    }
    return ComputePoolMeanAreaThreshold();
}

Return ModelPool::RecomputeAreaThresholds()
{
    m_area_threshold = ComputePoolAreaThreshold();
    if (m_area_threshold <= k_zero)
    {
        return Return::EmptyModel;
    }
    try
    {
        m_model_area_thresholds.assign(m_models.size(), m_area_threshold);
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }
    return Return::Success;
}

Return ModelPool::SetAppendedModelAreaThreshold(const std::size_t model_index)
{
    if (model_index >= m_models.size())
    {
        return Return::InvalidIndex;
    }
    if (m_model_area_thresholds.size() < m_models.size())
    {
        try
        {
            m_model_area_thresholds.resize(m_models.size(), m_area_threshold);
        }
        catch (const std::bad_alloc&)
        {
            return Return::OutOfMemory;
        }
    }
    m_model_area_thresholds[model_index] = m_area_threshold;
    return Return::Success;
}

Real ModelPool::ModelAreaThreshold(const std::size_t model_index) const
{
    if (model_index < m_model_area_thresholds.size() && m_model_area_thresholds[model_index] > k_zero)
    {
        return m_model_area_thresholds[model_index];
    }
    return m_area_threshold;
}

Return ModelPool::CountTriSubdivisions(const Vec3&  p1,
                                       const Vec3&  p2,
                                       const Vec3&  p3,
                                       const Real   area_threshold,
                                       std::size_t& subdivision_counter) const
{
    if (TriArea(p1, p2, p3) > area_threshold)
    {
        const Vec3 p12 = (p1 + p2) / 2;
        const Vec3 p23 = (p2 + p3) / 2;
        const Vec3 p31 = (p3 + p1) / 2;

        subdivision_counter++;

        return (CountTriSubdivisions(p1, p12, p31, area_threshold, subdivision_counter) == Return::Success &&
                CountTriSubdivisions(p12, p2, p23, area_threshold, subdivision_counter) == Return::Success &&
                CountTriSubdivisions(p31, p23, p3, area_threshold, subdivision_counter) == Return::Success &&
                CountTriSubdivisions(p12, p23, p31, area_threshold, subdivision_counter) == Return::Success)
                   ? Return::Success
                   : Return::OutOfMemory;
    }
    return Return::Success;
}

Return ModelPool::ComputeModelSubdivisionPlan(const build::Model& model,
                                              const Real          area_threshold,
                                              std::size_t&        model_subdivision_times,
                                              std::uint64_t&      subdivision_signature) const
{
    constexpr std::uint64_t fnv_offset = 14695981039346656037ull;
    constexpr std::uint64_t fnv_prime  = 1099511628211ull;
    auto                    hashValue  = [](std::uint64_t& hash, const auto& value)
    {
        const auto* bytes = reinterpret_cast<const unsigned char*>(&value);
        for (std::size_t i = 0; i < sizeof(value); ++i)
        {
            hash ^= static_cast<std::uint64_t>(bytes[i]);
            hash *= fnv_prime;
        }
    };

    model_subdivision_times     = 0;
    subdivision_signature       = fnv_offset;
    const std::size_t tri_count = model.OriginalTris().size();
    hashValue(subdivision_signature, tri_count);
    for (std::size_t tri_i = 0; tri_i < tri_count; ++tri_i)
    {
        const Tri&  tri                   = model.OriginalTris()[tri_i];
        std::size_t tri_subdivision_times = 0;
        Return      ret = CountTriSubdivisions(tri.p1, tri.p2, tri.p3, area_threshold, tri_subdivision_times);
        if (ret != Return::Success)
        {
            return ret;
        }
        model_subdivision_times += tri_subdivision_times;
        hashValue(subdivision_signature, tri_i);
        hashValue(subdivision_signature, tri_subdivision_times);
    }
    return Return::Success;
}

void ModelPool::DiscardTemporaryBuildData(build::Model& model) const
{
    std::vector<SubTri>().swap(model.m_tris);
    std::vector<std::size_t>().swap(model.m_leaf_subtri_indices);
}

// ============================================================================
// Build cache serialization
// ============================================================================

namespace
{
constexpr std::uint32_t       kBuildCacheVersion = 26;
constexpr std::array<char, 8> kBuildCacheMagic   = {'P', 'Q', 'P', 'C', 'C', '0', '8', '\0'};

struct BuildCacheHeader
{
    char          magic[8]{};
    std::uint32_t version                      = 0;
    std::uint32_t real_size                    = 0;
    std::uint32_t build_strategy               = 0;
    std::uint32_t reserved                     = 0;
    std::uint64_t source_content_hash          = 0;
    std::uint64_t reserved_source              = 0;
    std::uint64_t original_tri_count           = 0;
    std::uint64_t bv_count                     = 0;
    std::uint64_t leaf_count                   = 0;
    std::uint64_t subdivision_count            = 0;
    std::uint64_t subdivision_signature        = 0;
    std::uint64_t subtri_count                 = 0;
    std::uint64_t leaf_subtri_index_count      = 0;
    std::uint64_t b12v_rss_node_count          = 0;
    std::uint64_t m_b12v_rss_first_single_leaf = 0;
};

static_assert(std::is_trivially_copyable_v<BuildCacheHeader>);

template <typename T>
bool ReadExact(std::istream& is, T& value)
{
    return static_cast<bool>(is.read(reinterpret_cast<char*>(&value), static_cast<std::streamsize>(sizeof(T))));
}

template <typename T>
bool WriteExact(std::ostream& os, const T& value)
{
    return static_cast<bool>(os.write(reinterpret_cast<const char*>(&value), static_cast<std::streamsize>(sizeof(T))));
}

template <typename T>
bool ReadVector(std::istream& is, std::vector<T>& values)
{
    if (values.empty())
    {
        return true;
    }
    return static_cast<bool>(
        is.read(reinterpret_cast<char*>(values.data()), static_cast<std::streamsize>(values.size() * sizeof(T))));
}

template <typename T>
bool WriteVector(std::ostream& os, const std::vector<T>& values)
{
    if (values.empty())
    {
        return true;
    }
    return static_cast<bool>(os.write(reinterpret_cast<const char*>(values.data()),
                                      static_cast<std::streamsize>(values.size() * sizeof(T))));
}

bool ReadVec3(std::istream& is, Vec3& value)
{
    return ReadExact(is, value[0]) && ReadExact(is, value[1]) && ReadExact(is, value[2]);
}

bool WriteVec3(std::ostream& os, const Vec3& value)
{
    return WriteExact(os, value[0]) && WriteExact(os, value[1]) && WriteExact(os, value[2]);
}

bool ReadMat3(std::istream& is, Mat3& value)
{
    for (int r = 0; r < 3; ++r)
    {
        for (int c = 0; c < 3; ++c)
        {
            if (!ReadExact(is, value[r][c]))
            {
                return false;
            }
        }
    }
    return true;
}

bool WriteMat3(std::ostream& os, const Mat3& value)
{
    for (int r = 0; r < 3; ++r)
    {
        for (int c = 0; c < 3; ++c)
        {
            if (!WriteExact(os, value[r][c]))
            {
                return false;
            }
        }
    }
    return true;
}

bool ReadSubTri(std::istream& is, SubTri& value)
{
    return ReadVec3(is, value.p1) && ReadVec3(is, value.p2) && ReadVec3(is, value.p3) &&
           ReadExact(is, value.original_tri_id);
}

bool WriteSubTri(std::ostream& os, const SubTri& value)
{
    return WriteVec3(os, value.p1) && WriteVec3(os, value.p2) && WriteVec3(os, value.p3) &&
           WriteExact(os, value.original_tri_id);
}

bool WriteBV(std::ostream& os, const build::BV& value)
{
    const std::array<Real, 2>& length      = value.L();
    const Real                 radius      = value.Radius();
    const Real                 size        = value.Size();
    const int                  first_child = value.FirstChildRaw();
    const int                  tri0        = value.Tri0();
    const int                  tri1        = value.Tri1();
    const bool shape_written = WriteMat3(os, value.R()) && WriteVec3(os, value.Tr()) && WriteExact(os, length[0]) &&
                               WriteExact(os, length[1]) && WriteExact(os, radius) && WriteExact(os, size);
    return shape_written && WriteExact(os, first_child) && WriteExact(os, tri0) && WriteExact(os, tri1);
}

bool ValidateB12vTriIds(const std::vector<build::BV>& bvs,
                        const std::size_t             original_tri_count,
                        const std::size_t             first_single_leaf)
{
    if (bvs.empty() || bvs.size() % 2 == 0)
    {
        return false;
    }

    const std::size_t first_leaf = bvs.size() / 2;
    if (first_single_leaf < first_leaf || first_single_leaf > bvs.size())
    {
        return false;
    }
    for (std::size_t node = 0; node < first_leaf; ++node)
    {
        if (bvs[node].FirstChildRaw() != -1 || bvs[node].Tri0() != -1 || bvs[node].Tri1() != -1)
        {
            return false;
        }
    }

    for (std::size_t node = first_leaf; node < bvs.size(); ++node)
    {
        const int tri0 = bvs[node].Tri0();
        const int tri1 = bvs[node].Tri1();
        if (bvs[node].FirstChildRaw() != -1 || tri0 < 0 || static_cast<std::size_t>(tri0) >= original_tri_count ||
            tri1 < -1 || (tri1 >= 0 && static_cast<std::size_t>(tri1) >= original_tri_count))
        {
            return false;
        }
    }
    return true;
}

std::filesystem::path NormalizePath(const std::filesystem::path& path)
{
    std::error_code       ec;
    std::filesystem::path normalized = std::filesystem::weakly_canonical(path, ec);
    if (ec)
    {
        ec.clear();
        normalized = std::filesystem::absolute(path, ec);
        if (ec)
        {
            return path.lexically_normal();
        }
    }
    return normalized.lexically_normal();
}

struct SourceFingerprint
{
    std::filesystem::path normalized_path;
    std::uint64_t         content_hash = 0;
};

bool MakeSourceFingerprint(const std::filesystem::path& filename, SourceFingerprint& fingerprint)
{
    fingerprint.normalized_path = NormalizePath(filename);

    std::ifstream file(fingerprint.normalized_path, std::ios::binary);
    if (!file)
    {
        return false;
    }

    constexpr std::uint64_t     fnv_offset = 14695981039346656037ull;
    constexpr std::uint64_t     fnv_prime  = 1099511628211ull;
    std::uint64_t               hash       = fnv_offset;
    std::array<char, 64 * 1024> buffer{};
    while (file)
    {
        file.read(buffer.data(), static_cast<std::streamsize>(buffer.size()));
        const std::streamsize bytes_read = file.gcount();
        for (std::streamsize i = 0; i < bytes_read; ++i)
        {
            hash ^= static_cast<std::uint64_t>(static_cast<unsigned char>(buffer[static_cast<std::size_t>(i)]));
            hash *= fnv_prime;
        }
    }
    fingerprint.content_hash = hash;
    return true;
}
} // namespace

bool ModelPool::TryLoadBuildCache(build::Model&          model,
                                  const ModelSourceInfo& source,
                                  const std::size_t      model_subdivision_times,
                                  const std::uint64_t    subdivision_signature)
{
    if (!source.from_file)
    {
        return false;
    }

    SourceFingerprint fingerprint;
    if (!MakeSourceFingerprint(source.filename, fingerprint))
    {
        return false;
    }

    std::filesystem::path cache_path;
    if (!pqss::build_cache_path::MakeBuildCacheFilePath(fingerprint.content_hash,
                                                        model.m_build_strategy,
                                                        std::is_same_v<Real, float>,
                                                        m_build_cache_tag,
                                                        model_subdivision_times,
                                                        subdivision_signature,
                                                        cache_path))
    {
        return false;
    }

    std::ifstream ifs(cache_path, std::ios::binary);
    if (!ifs)
    {
        return false;
    }

    BuildCacheHeader header;
    if (!ReadExact(ifs, header))
    {
        return false;
    }
    if (!std::equal(kBuildCacheMagic.begin(), kBuildCacheMagic.end(), std::begin(header.magic)))
    {
        return false;
    }
    if (header.version != kBuildCacheVersion || header.real_size != sizeof(Real) ||
        header.build_strategy != static_cast<std::uint32_t>(model.m_build_strategy) ||
        header.source_content_hash != fingerprint.content_hash ||
        header.original_tri_count != model.m_original_tris.size() ||
        header.subdivision_count != model_subdivision_times || header.subdivision_signature != subdivision_signature)
    {
        return false;
    }

    try
    {
        std::vector<build::BV>   cached_b(static_cast<std::size_t>(header.bv_count));
        std::vector<SubTri>      cached_tris(static_cast<std::size_t>(header.subtri_count));
        std::vector<std::size_t> cached_leaf_subtri_indices(static_cast<std::size_t>(header.leaf_subtri_index_count));
        for (build::BV& bv : cached_b)
        {
            if (!ReadMat3(ifs, bv.m_R) || !ReadVec3(ifs, bv.m_tr) || !ReadExact(ifs, bv.m_l[0]) ||
                !ReadExact(ifs, bv.m_l[1]) || !ReadExact(ifs, bv.m_r) || !ReadExact(ifs, bv.m_size) ||
                !ReadExact(ifs, bv.m_first_child) || !ReadExact(ifs, bv.m_tri0) || !ReadExact(ifs, bv.m_tri1))
            {
                return false;
            }
        }
        for (SubTri& tri : cached_tris)
        {
            if (!ReadSubTri(ifs, tri))
            {
                return false;
            }
        }
        if (!ReadVector(ifs, cached_leaf_subtri_indices))
        {
            return false;
        }
        if (header.leaf_count != (header.bv_count == 0 ? 0 : (header.bv_count + 1) / 2))
        {
            return false;
        }
        model.m_b                   = std::move(cached_b);
        model.m_tris                = std::move(cached_tris);
        model.m_leaf_subtri_indices = std::move(cached_leaf_subtri_indices);
        model.m_next_bv_idx         = model.m_b.size();
        if (m_build_b12v_enabled && header.b12v_rss_node_count > 0)
        {
            std::vector<build::BV> cached_b12v(static_cast<std::size_t>(header.b12v_rss_node_count));
            for (build::BV& bv : cached_b12v)
            {
                if (!ReadMat3(ifs, bv.m_R) || !ReadVec3(ifs, bv.m_tr) || !ReadExact(ifs, bv.m_l[0]) ||
                    !ReadExact(ifs, bv.m_l[1]) || !ReadExact(ifs, bv.m_r) || !ReadExact(ifs, bv.m_size) ||
                    !ReadExact(ifs, bv.m_first_child) || !ReadExact(ifs, bv.m_tri0) || !ReadExact(ifs, bv.m_tri1))
                {
                    return false;
                }
            }
            if (!ValidateB12vTriIds(cached_b12v,
                                    static_cast<std::size_t>(header.original_tri_count),
                                    static_cast<std::size_t>(header.m_b12v_rss_first_single_leaf)))
            {
                return false;
            }
            model.m_b12v_rss_first_single_leaf = static_cast<std::size_t>(header.m_b12v_rss_first_single_leaf);
            model.m_b12v_rss_b                 = std::move(cached_b12v);
        }
        return true;
    }
    catch (const std::bad_alloc&)
    {
        return false;
    }
}

bool ModelPool::SaveBuildCache(const build::Model&    model,
                               const ModelSourceInfo& source,
                               const std::size_t      model_subdivision_times,
                               const std::uint64_t    subdivision_signature) const
{
    if (!source.from_file)
    {
        return false;
    }

    SourceFingerprint fingerprint;
    if (!MakeSourceFingerprint(source.filename, fingerprint))
    {
        return false;
    }

    std::filesystem::path cache_path;
    if (!pqss::build_cache_path::MakeBuildCacheFilePath(fingerprint.content_hash,
                                                        model.m_build_strategy,
                                                        std::is_same_v<Real, float>,
                                                        m_build_cache_tag,
                                                        model_subdivision_times,
                                                        subdivision_signature,
                                                        cache_path))
    {
        return false;
    }

    std::error_code ec;
    std::filesystem::create_directories(cache_path.parent_path(), ec);
    if (ec)
    {
        return false;
    }

    std::ofstream ofs(cache_path, std::ios::binary | std::ios::trunc);
    if (!ofs)
    {
        return false;
    }

    BuildCacheHeader header;
    std::copy(kBuildCacheMagic.begin(), kBuildCacheMagic.end(), std::begin(header.magic));
    header.version                 = kBuildCacheVersion;
    header.real_size               = sizeof(Real);
    header.build_strategy          = static_cast<std::uint32_t>(model.m_build_strategy);
    header.source_content_hash     = fingerprint.content_hash;
    header.original_tri_count      = static_cast<std::uint64_t>(model.m_original_tris.size());
    header.bv_count                = static_cast<std::uint64_t>(model.m_b.size());
    header.leaf_count              = static_cast<std::uint64_t>(model.m_b.empty() ? 0 : (model.m_b.size() + 1) / 2);
    header.subdivision_count       = static_cast<std::uint64_t>(model_subdivision_times);
    header.subdivision_signature   = subdivision_signature;
    header.subtri_count            = static_cast<std::uint64_t>(model.m_tris.size());
    header.leaf_subtri_index_count = static_cast<std::uint64_t>(model.m_leaf_subtri_indices.size());
    header.b12v_rss_node_count     = static_cast<std::uint64_t>(model.m_b12v_rss_b.size());
    header.m_b12v_rss_first_single_leaf = static_cast<std::uint64_t>(model.m_b12v_rss_first_single_leaf);

    if (!WriteExact(ofs, header))
    {
        return false;
    }
    for (const build::BV& bv : model.m_b)
    {
        if (!WriteBV(ofs, bv))
        {
            return false;
        }
    }
    for (const SubTri& tri : model.m_tris)
    {
        if (!WriteSubTri(ofs, tri))
        {
            return false;
        }
    }
    if (!WriteVector(ofs, model.m_leaf_subtri_indices))
    {
        return false;
    }
    for (const build::BV& bv : model.m_b12v_rss_b)
    {
        if (!WriteBV(ofs, bv))
        {
            return false;
        }
    }
    return true;
}

// ============================================================================
// Build model data helpers
// ============================================================================

Return ModelPool::BuildModelDataWithCurrentThreshold(build::Model& model,
                                                     const Real    area_threshold,
                                                     std::size_t&  model_subdivision_times) const
{
    model_subdivision_times = 0;
    try
    {
        model.m_tris.clear();
        model.m_b.clear();
        model.m_leaf_subtri_indices.clear();
        model.m_tris.reserve(model.m_original_tris.size());
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    for (std::size_t tri_i = 0; tri_i < model.m_original_tris.size(); ++tri_i)
    {
        const Tri& original_tri = model.m_original_tris[tri_i];
        Return     ret          = AddTriRecurseCounting(model,
                                                        original_tri.p1,
                                                        original_tri.p2,
                                                        original_tri.p3,
                                                        tri_i,
                                                        area_threshold,
                                                        model_subdivision_times);
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    return model.EndPool(true);
}

Return ModelPool::BuildModelDataWithoutSubdivision(build::Model& model, std::size_t& model_subdivision_times) const
{
    model_subdivision_times = 0;
    try
    {
        model.m_tris.clear();
        model.m_b.clear();
        model.m_leaf_subtri_indices.clear();
        model.m_tris.reserve(model.m_original_tris.size());
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    for (std::size_t tri_i = 0; tri_i < model.m_original_tris.size(); ++tri_i)
    {
        const Tri& original_tri = model.m_original_tris[tri_i];
        Return     ret          = model.AddTri(original_tri.p1, original_tri.p2, original_tri.p3, tri_i);
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    return model.EndPool(true);
}

Return ModelPool::BuildModelWithCurrentThreshold(build::Model&          model,
                                                 const ModelSourceInfo& source,
                                                 const Real             area_threshold,
                                                 std::size_t&           model_subdivision_times)
{
    std::uint64_t subdivision_signature = 0;
    Return        ret                   = Return::Success;
    if (m_subdivision_enabled)
    {
        ret = ComputeModelSubdivisionPlan(model, area_threshold, model_subdivision_times, subdivision_signature);
        if (ret != Return::Success)
        {
            return ret;
        }
    }
    else
    {
        ret = ComputeModelSubdivisionPlan(model,
                                          std::numeric_limits<Real>::max(),
                                          model_subdivision_times,
                                          subdivision_signature);
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    if (source.from_file)
    {
        m_build_cache_stats.eligible_models++;
        Timer phase_timer;
        phase_timer.Start();
        const bool cache_hit = TryLoadBuildCache(model, source, model_subdivision_times, subdivision_signature);
        phase_timer.Pause();
        m_build_cache_stats.cache_load_ms += phase_timer.ElapsedMs();
        if (cache_hit)
        {
            m_build_cache_stats.cache_hits++;
            if (m_build_b12v_enabled && model.m_b12v_rss_b.empty())
            {
                ret = model.BuildB12vRss();
                if (ret != Return::Success)
                {
                    return ret;
                }
            }
            DiscardTemporaryBuildData(model);
            return Return::Success;
        }
        m_build_cache_stats.cache_misses++;
    }

    Timer build_timer;
    build_timer.Start();
    if (m_subdivision_enabled)
    {
        ret = BuildModelDataWithCurrentThreshold(model, area_threshold, model_subdivision_times);
    }
    else
    {
        ret = BuildModelDataWithoutSubdivision(model, model_subdivision_times);
    }
    build_timer.Pause();
    m_build_cache_stats.cache_build_ms += build_timer.ElapsedMs();
    if (ret != Return::Success)
    {
        return ret;
    }

    if (m_build_b12v_enabled)
    {
        ret = model.BuildB12vRss();
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    if (source.from_file)
    {
        Timer phase_timer;
        phase_timer.Start();
        SaveBuildCache(model, source, model_subdivision_times, subdivision_signature);
        phase_timer.Pause();
        m_build_cache_stats.cache_save_ms += phase_timer.ElapsedMs();
    }
    DiscardTemporaryBuildData(model);
    return Return::Success;
}

Return ModelPool::BuildAppendedModel(const std::size_t model_index)
{
    if (model_index >= m_models.size())
    {
        return Return::InvalidIndex;
    }
    Return threshold_ret = SetAppendedModelAreaThreshold(model_index);
    if (threshold_ret != Return::Success)
    {
        return threshold_ret;
    }

    std::size_t model_subdivision_times = 0;
    Return      ret                     = BuildModelWithCurrentThreshold(m_models[model_index],
                                                                         m_model_sources[model_index],
                                                                         ModelAreaThreshold(model_index),
                                                                         model_subdivision_times);
    if (ret == Return::Success)
    {
        m_subdivision_times += model_subdivision_times;
    }
    return ret;
}

Return ModelPool::AddTriRecurse(build::Model&     model,
                                const Vec3&       p1,
                                const Vec3&       p2,
                                const Vec3&       p3,
                                const std::size_t original_tri_id)
{
    std::size_t subdivision_count = 0;
    Return      ret = AddTriRecurseCounting(model, p1, p2, p3, original_tri_id, m_area_threshold, subdivision_count);
    if (ret == Return::Success)
    {
        m_subdivision_times += subdivision_count;
    }
    return ret;
}

Return ModelPool::AddTriRecurseCounting(build::Model&     model,
                                        const Vec3&       p1,
                                        const Vec3&       p2,
                                        const Vec3&       p3,
                                        const std::size_t original_tri_id,
                                        const Real        area_threshold,
                                        std::size_t&      subdivision_counter) const
{
    if (TriArea(p1, p2, p3) > area_threshold)
    {
        const Vec3 p12 = (p1 + p2) / 2;
        const Vec3 p23 = (p2 + p3) / 2;
        const Vec3 p31 = (p3 + p1) / 2;

        subdivision_counter++;

        return (AddTriRecurseCounting(model, p1, p12, p31, original_tri_id, area_threshold, subdivision_counter) ==
                    Return::Success &&
                AddTriRecurseCounting(model, p12, p2, p23, original_tri_id, area_threshold, subdivision_counter) ==
                    Return::Success &&
                AddTriRecurseCounting(model, p31, p23, p3, original_tri_id, area_threshold, subdivision_counter) ==
                    Return::Success &&
                AddTriRecurseCounting(model, p12, p23, p31, original_tri_id, area_threshold, subdivision_counter) ==
                    Return::Success)
                   ? Return::Success
                   : Return::OutOfMemory;
    }

    return model.AddTri(p1, p2, p3, original_tri_id);
}

// ============================================================================
// Build API: LoadModelFromFile
// ============================================================================

Return ModelPool::LoadModelFromFile(const std::string_view filename, const BuildStrategy build_strategy)
{
    const bool append_to_processed_pool = m_build_state == PoolBuildState::Processed;
    const std::size_t original_pool_size = m_models.size();

    build::Model new_model;
    Return       ret = new_model.LoadModelFromFile(filename, build_strategy);
    if (ret != Return::Success)
    {
        std::cout << "Warning: failed to read model file; the model pool is unchanged.\n";
        return ret;
    }

    try
    {
        m_models.emplace_back(std::move(new_model));
        m_model_sources.push_back(ModelSourceInfo{true, std::filesystem::u8path(filename)});
        m_model_area_thresholds.push_back(k_zero);
        m_rss_leaf_collision_models.push_back(0);
    }
    catch (const std::bad_alloc&)
    {
        m_models.resize(original_pool_size);
        m_model_sources.resize(original_pool_size);
        m_model_area_thresholds.resize(original_pool_size);
        m_rss_leaf_collision_models.resize(original_pool_size);
        return Return::OutOfMemory;
    }

    if (append_to_processed_pool)
    {
        ret = BuildAppendedModel(m_models.size() - 1);
        if (ret != Return::Success)
        {
            m_models.pop_back();
            m_model_sources.pop_back();
            m_model_area_thresholds.pop_back();
            m_rss_leaf_collision_models.pop_back();
            return ret;
        }
        m_build_state = PoolBuildState::Processed;
        return Return::Success;
    }

    m_build_state = PoolBuildState::Ready;
    return Return::Success;
}

Return ModelPool::Reserve(const std::size_t num_models)
{
    try
    {
        m_models.reserve(num_models);
        m_model_sources.reserve(num_models);
        m_model_area_thresholds.reserve(num_models);
        m_rss_leaf_collision_models.reserve(num_models);
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    if (m_build_state != PoolBuildState::Processed)
    {
        m_build_state = PoolBuildState::Ready;
    }
    return Return::Success;
}

Return ModelPool::ReserveTris(const std::size_t num_tris)
{
    if (m_build_state != PoolBuildState::Building)
    {
        return Return::NotBegun;
    }

    return m_models.back().Reserve(num_tris);
}

Return ModelPool::NewModel(const std::size_t num_tris, const BuildStrategy build_strategy)
{
    if (num_tris == 0)
    {
        return Return::InvalidIndex;
    }
    if (m_build_state == PoolBuildState::Empty)
    {
        std::cout << "Warning: no capacity was reserved; reserving space for 8 models.\n";
        Return ret = Reserve(8);
        if (ret != Return::Success)
        {
            return ret;
        }
    }
    if (m_build_state == PoolBuildState::Building)
    {
        std::cout << "Warning: finalizing the current model automatically.\n";
        Return ret = EndModel();
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    const bool append_to_processed_pool = m_build_state == PoolBuildState::Processed;
    const std::size_t original_pool_size = m_models.size();

    try
    {
        m_models.emplace_back();
        m_model_sources.emplace_back();
        m_model_area_thresholds.push_back(k_zero);
        m_rss_leaf_collision_models.push_back(0);
    }
    catch (const std::bad_alloc&)
    {
        m_models.resize(original_pool_size);
        m_model_sources.resize(original_pool_size);
        m_model_area_thresholds.resize(original_pool_size);
        m_rss_leaf_collision_models.resize(original_pool_size);
        return Return::OutOfMemory;
    }

    m_build_state                    = PoolBuildState::Building;
    m_append_with_current_threshold  = append_to_processed_pool;
    m_models.back().m_build_strategy = build_strategy;
    return Return::Success;
}

Return ModelPool::AddTri(const Vec3& p1, const Vec3& p2, const Vec3& p3)
{
    if (m_build_state == PoolBuildState::Ready || m_build_state == PoolBuildState::Empty ||
        m_build_state == PoolBuildState::Processed)
    {
        std::cout << "Warning: no model is being built; creating one automatically.\n";
        Return ret = NewModel(8, BuildStrategy::Optimized);
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    return m_models.back().AddOriginalTri(p1, p2, p3);
}

Return ModelPool::ResetBuildStrategy(const BuildStrategy build_strategy)
{
    if (m_build_state != PoolBuildState::Building)
    {
        return Return::NotBegun;
    }

    m_models.back().m_build_strategy = build_strategy;
    return Return::Success;
}

Return ModelPool::EndModel()
{
    if (m_build_state != PoolBuildState::Building)
    {
        return Return::NotBegun;
    }

    Return ret = m_models.back().EndModel();
    if (ret != Return::Success)
    {
        return ret;
    }
    if (m_append_with_current_threshold)
    {
        ret = BuildAppendedModel(m_models.size() - 1);
        if (ret != Return::Success)
        {
            return ret;
        }
        m_append_with_current_threshold = false;
        m_build_state                   = PoolBuildState::Processed;
        return Return::Success;
    }

    m_build_state = PoolBuildState::Ready;
    return Return::Success;
}

Return ModelPool::DeleteModel(const std::size_t model_id)
{
    if (m_build_state == PoolBuildState::Processed)
    {
        return Return::AlreadyBuilt;
    }
    if (model_id >= m_models.size())
    {
        return Return::InvalidIndex;
    }
    if (m_build_state == PoolBuildState::Building && model_id == m_models.size() - 1)
    {
        std::cout << "Warning: finalizing the current model before deletion.\n";
        Return ret = EndModel();
        if (ret != Return::Success)
        {
            return ret;
        }
    }

    m_models.erase(m_models.begin() + static_cast<std::ptrdiff_t>(model_id));
    m_model_sources.erase(m_model_sources.begin() + static_cast<std::ptrdiff_t>(model_id));
    m_model_area_thresholds.erase(m_model_area_thresholds.begin() + static_cast<std::ptrdiff_t>(model_id));
    m_rss_leaf_collision_models.erase(
        m_rss_leaf_collision_models.begin() + static_cast<std::ptrdiff_t>(model_id));
    return Return::Success;
}

Return ModelPool::EndPool()
{
    if (m_build_state == PoolBuildState::Processed)
    {
        return Return::Success;
    }
    if (m_build_state == PoolBuildState::Empty)
    {
        return Return::EmptyModel;
    }
    if (m_build_state == PoolBuildState::Building)
    {
        std::cout << "Warning: finalizing the current model before building the pool.\n";
        Return ret = EndModel();
        if (ret != Return::Success)
        {
            return ret;
        }
        if (m_build_state == PoolBuildState::Processed)
        {
            return Return::Success;
        }
    }

    m_build_cache_stats = {};
    m_subdivision_times = 0;

    Timer phase_timer;
    phase_timer.Start();

    Return threshold_ret = RecomputeAreaThresholds();
    if (threshold_ret != Return::Success)
    {
        return threshold_ret;
    }

    phase_timer.Pause();
    m_build_cache_stats.area_threshold_ms = phase_timer.ElapsedMs();

    for (std::size_t model_i = 0; model_i < m_models.size(); ++model_i)
    {
        build::Model&          model                   = m_models[model_i];
        const ModelSourceInfo& source                  = m_model_sources[model_i];
        std::size_t            model_subdivision_times = 0;

        Return ret =
            BuildModelWithCurrentThreshold(model, source, ModelAreaThreshold(model_i), model_subdivision_times);
        if (ret != Return::Success)
        {
            return ret;
        }
        m_subdivision_times += model_subdivision_times;
    }
    m_append_with_current_threshold = false;
    m_build_state                   = PoolBuildState::Processed;
    return Return::Success;
}

Return ModelPool::Clear()
{
    m_build_state                   = PoolBuildState::Empty;
    m_area_threshold                = k_zero;
    m_area_threshold_override       = k_zero;
    m_subdivision_times             = 0;
    m_append_with_current_threshold = false;
    m_build_cache_stats             = {};

    std::vector<build::Model>().swap(m_models);
    std::vector<ModelSourceInfo>().swap(m_model_sources);
    std::vector<Real>().swap(m_model_area_thresholds);
    std::vector<std::uint8_t>().swap(m_rss_leaf_collision_models);

    return Return::Success;
}

Return ModelPool::SetSubdivisionEnabled(const bool enabled)
{
    if (m_build_state == PoolBuildState::Processed)
    {
        return Return::AlreadyBuilt;
    }
    m_subdivision_enabled = enabled;
    return Return::Success;
}

Return ModelPool::SetAreaThresholdOverride(const Real area_threshold)
{
    if (m_build_state == PoolBuildState::Processed ||
        m_build_state == PoolBuildState::Building)
    {
        return Return::AlreadyBuilt;
    }
    if (area_threshold < k_zero)
    {
        return Return::UnknownError;
    }
    m_area_threshold_override = area_threshold;
    return Return::Success;
}

Return ModelPool::SetB12vBuildEnabled(const bool enabled)
{
    if (m_build_state == PoolBuildState::Processed)
    {
        return Return::AlreadyBuilt;
    }
    m_build_b12v_enabled = enabled;
    return Return::Success;
}

Return ModelPool::SetBuildCacheTag(const std::string_view cache_tag)
{
    if (m_build_state == PoolBuildState::Processed)
    {
        return Return::AlreadyBuilt;
    }
    if (cache_tag.empty())
    {
        return Return::InvalidIndex;
    }
    try
    {
        m_build_cache_tag = cache_tag;
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }
    return Return::Success;
}

Return ModelPool::ReplaceModelBvhWithRssLeaves(const std::size_t model_id,
                                               const std::vector<RssBvhNode>& nodes)
{
    if (m_build_state != PoolBuildState::Processed)
        return Return::NotProccessed;
    if (model_id >= m_models.size() || nodes.empty())
        return Return::InvalidIndex;

    std::vector<build::BV> packed;
    std::vector<std::uint8_t> visited;
    try
    {
        packed.resize(nodes.size());
        visited.resize(nodes.size(), 0);
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    std::size_t next_node = 1;
    std::function<bool(std::size_t, std::size_t, const Mat3&, const Vec3&)> pack;
    pack = [&](const std::size_t source_index,
               const std::size_t destination_index,
               const Mat3& parent_R,
               const Vec3& parent_T)
    {
        if (source_index >= nodes.size() || destination_index >= packed.size() || visited[source_index] != 0)
            return false;
        visited[source_index] = 1;
        const RssBvhNode& source = nodes[source_index];
        if (!std::isfinite(source.radius) || source.radius < k_zero ||
            !std::isfinite(source.lengths[0]) || source.lengths[0] < k_zero ||
            !std::isfinite(source.lengths[1]) || source.lengths[1] < k_zero)
            return false;

        build::BV& destination = packed[destination_index];
        destination.m_R = MTxM(parent_R, source.R);
        destination.m_tr = MTxV(parent_R, source.T - parent_T);
        destination.m_l = source.lengths;
        destination.m_r = source.radius;
        destination.m_size = source.size;
        if (source.leaf)
        {
            destination.SetLeaf(0);
            return true;
        }
        if (source.left >= nodes.size() || source.right >= nodes.size() ||
            source.left == source.right || next_node + 1 >= packed.size())
            return false;
        const std::size_t first_child = next_node;
        next_node += 2;
        destination.SetInternal(first_child);
        return pack(source.left, first_child, source.R, source.T) &&
               pack(source.right, first_child + 1, source.R, source.T);
    };

    if (!pack(0, 0, Meident(), Veident()) || next_node != packed.size() ||
        std::find(visited.begin(), visited.end(), 0) != visited.end())
        return Return::InvalidIndex;

    build::Model& model = m_models[model_id];
    model.m_b = std::move(packed);
    model.m_b12v_rss_b.clear();
    model.m_b12v_rss_first_single_leaf = 0;
    model.m_tris.clear();
    model.m_leaf_subtri_indices.clear();
    m_rss_leaf_collision_models[model_id] = 1;
    return Return::Success;
}

// ============================================================================
// Query API
// ============================================================================

bool ModelPool::Processed() const
{
    return m_build_state == PoolBuildState::Processed;
}

std::size_t ModelPool::PoolSize() const
{
    if (m_build_state != PoolBuildState::Processed)
    {
        std::cout << "Warning: model count is unavailable before the pool build completes.\n";
        return 0;
    }
    return m_models.size();
}

std::size_t ModelPool::Info(const bool print) const
{
    if (m_build_state != PoolBuildState::Processed)
    {
        std::cout << "Warning: model pool information is unavailable before the build completes.\n";
        return 0;
    }

    std::size_t total_mem_usage = sizeof(ModelPool);
    for (const build::Model& model : m_models)
    {
        total_mem_usage += model.MemUsage();
    }

    if (print)
    {
        std::cout << std::format(
            "\n---------- Model pool information ----------\nAddress: {}\tMemory: {} bytes\tModels: {}\n",
            static_cast<const void*>(this),
            total_mem_usage,
            m_models.size());

        for (std::size_t i = 0; i < m_models.size(); ++i)
        {
            const build::Model& model = m_models[i];
            std::cout << std::format(
                "Model {}: address: {}\tmemory: {} bytes\toriginal triangles: {}\tbuilt triangles: {}\n",
                i,
                static_cast<const void*>(&model),
                model.MemUsage(),
                model.m_original_tris.size(),
                model.m_b.empty() ? 0 : (model.m_b.size() + 1) / 2);
        }
        std::cout << "------------------------------\n\n";
    }

    return total_mem_usage;
}

Real ModelPool::MaxArea() const
{
    if (m_build_state != PoolBuildState::Processed)
    {
        std::cout << "Warning: area threshold is unavailable before the pool build completes.\n";
        return k_zero;
    }
    return m_area_threshold;
}

std::size_t ModelPool::SubdivisionTimes() const
{
    if (m_build_state != PoolBuildState::Processed)
    {
        std::cout << "Warning: subdivision count is unavailable before the pool build completes.\n";
        return 0;
    }
    return m_subdivision_times;
}

bool ModelPool::SubdivisionEnabled() const
{
    return m_subdivision_enabled;
}

bool ModelPool::B12vBuildEnabled() const
{
    return m_build_b12v_enabled;
}

std::string_view ModelPool::BuildCacheTag() const
{
    return m_build_cache_tag;
}

const ModelPool::BuildCacheStats& ModelPool::GetBuildCacheStats() const
{
    return m_build_cache_stats;
}

std::vector<Real> ModelPool::ModelAreaThresholds() const
{
    return m_model_area_thresholds;
}

std::vector<ModelPool::ModelStats> ModelPool::GetModelStats() const
{
    std::vector<ModelStats> stats;
    stats.reserve(m_models.size());
    for (const build::Model& model : m_models)
    {
        const std::size_t original_tri_count = model.OriginalTris().size();
        const std::size_t built_tri_count    = model.Bvs().empty() ? 0 : (model.Bvs().size() + 1) / 2;
        const std::size_t subdivision_count =
            built_tri_count >= original_tri_count ? (built_tri_count - original_tri_count) / 3 : 0;
        stats.push_back(ModelStats{original_tri_count, built_tri_count, subdivision_count, model.Bvs().size()});
    }
    return stats;
}

std::vector<ModelPool::BvhSignature> ModelPool::BvhSignatures() const
{
    constexpr std::uint64_t fnv_offset = 14695981039346656037ull;
    constexpr std::uint64_t fnv_prime  = 1099511628211ull;

    auto hashBytes = [](std::uint64_t& h, const void* data, const std::size_t size)
    {
        const auto* bytes = static_cast<const unsigned char*>(data);
        for (std::size_t i = 0; i < size; ++i)
        {
            h ^= static_cast<std::uint64_t>(bytes[i]);
            h *= fnv_prime;
        }
    };
    auto hashValue = [&hashBytes](std::uint64_t& h, const auto& value)
    {
        hashBytes(h, &value, sizeof(value));
    };

    std::vector<BvhSignature> signatures;
    signatures.reserve(m_models.size());
    for (const build::Model& model : m_models)
    {
        BvhSignature                  signature;
        const std::vector<build::BV>& bvs = model.Bvs();
        signature.original_tri_count      = model.OriginalTris().size();
        signature.bv_count                = bvs.size();
        signature.built_tri_count         = bvs.empty() ? 0 : (bvs.size() + 1) / 2;
        signature.subdivision_count       = signature.built_tri_count >= signature.original_tri_count
                                                ? (signature.built_tri_count - signature.original_tri_count) / 3
                                                : 0;
        signature.topology_signature      = fnv_offset;
        signature.geometry_signature      = fnv_offset;
        signature.leaf_signature          = fnv_offset;

        hashValue(signature.topology_signature, signature.bv_count);
        hashValue(signature.geometry_signature, signature.bv_count);
        hashValue(signature.leaf_signature, signature.bv_count);

        for (std::size_t i = 0; i < bvs.size(); ++i)
        {
            const build::BV&    bv        = bvs[i];
            const bool          is_leaf   = bv.IsLeaf();
            const unsigned char node_kind = is_leaf ? 1 : 0;
            hashValue(signature.topology_signature, i);
            hashValue(signature.topology_signature, node_kind);
            if (is_leaf)
            {
                ++signature.leaf_count;
                signature.leaf_size_sum += bv.Size();
                const std::size_t tri_index = bv.TriIndex();
                hashValue(signature.leaf_signature, i);
                hashValue(signature.leaf_signature, tri_index);
            }
            else
            {
                ++signature.internal_count;
                signature.internal_size_sum += bv.Size();
                const std::size_t first_child = bv.FirstChild();
                hashValue(signature.topology_signature, first_child);
            }

            signature.total_size += bv.Size();
            hashValue(signature.geometry_signature, i);
            hashValue(signature.geometry_signature, bv.R());
            hashValue(signature.geometry_signature, bv.Tr());
            hashValue(signature.geometry_signature, bv.L());
            const Real radius = bv.Radius();
            const Real size   = bv.Size();
            hashValue(signature.geometry_signature, radius);
            hashValue(signature.geometry_signature, size);
        }

        if (!bvs.empty())
        {
            signature.root_size = bvs.front().Size();

            struct StackItem
            {
                std::size_t node  = 0;
                std::size_t depth = 0;
            };

            std::vector<StackItem> stack;
            stack.push_back(StackItem{0, 1});
            while (!stack.empty())
            {
                const StackItem item = stack.back();
                stack.pop_back();
                if (item.node >= bvs.size())
                {
                    continue;
                }
                signature.max_depth = std::max(signature.max_depth, item.depth);
                const build::BV& bv = bvs[item.node];
                if (!bv.IsLeaf())
                {
                    const std::size_t first_child = bv.FirstChild();
                    stack.push_back(StackItem{first_child, item.depth + 1});
                    stack.push_back(StackItem{first_child + 1, item.depth + 1});
                }
            }
        }

        signatures.push_back(signature);
    }
    return signatures;
}

const std::vector<build::Model>& ModelPool::Models() const
{
    return m_models;
}

// ============================================================================
// Collision detection
// ============================================================================

Return ModelPool::Query(QueryResult&      res,
                        const std::size_t model_id_1,
                        const Mat3&       R1,
                        const Vec3&       T1,
                        const std::size_t model_id_2,
                        const Mat3&       R2,
                        const Vec3&       T2,
                        const Real        tolerance) const
{
    Timer timer;
    timer.Start();

    if (m_build_state != PoolBuildState::Processed)
    {
        return Return::NotProccessed;
    }
    if (model_id_1 >= m_models.size() || model_id_2 >= m_models.size())
    {
        return Return::InvalidIndex;
    }

    res.m_R = MTxM(R1, R2);
    res.m_T = MTxVmV(R1, T2, T1);

    res.m_tolerance             = tolerance;
    res.m_num_bv_tests          = 0;
    res.m_num_tri_tests         = 0;
    res.m_closer_than_tolerance = 0;

    const build::Model& o1 = m_models[model_id_1];
    const build::Model& o2 = m_models[model_id_2];

    const build::BV& bv1 = o1.Bvs()[0];
    const build::BV& bv2 = o2.Bvs()[0];

    Mat3 R = MTxM(bv1.R(), MxM(res.R(), bv2.R()));
    Vec3 T = MTxV(bv1.R(), MxVpV(res.R(), bv2.Tr(), res.T()) - bv1.Tr());

    Real d                = detail::RectDistSq(R, T, bv1.L(), bv2.L());
    Real tolerance_plus_r = res.Tolerance() + bv1.Radius() + bv2.Radius();
    if (d <= tolerance_plus_r * tolerance_plus_r)
    {
        QueryRecurse(res, R, T, o1, 0, o2, 0,
                     m_rss_leaf_collision_models[model_id_1] != 0 ||
                         m_rss_leaf_collision_models[model_id_2] != 0);
    }

    res.m_p2            = MTxVmV(res.R(), res.P2().value_or(Vec3()), res.T());
    res.m_query_time_ms = timer.ElapsedMs();

    return Return::Success;
}

void ModelPool::QueryRecurse(QueryResult&        res,
                             const Mat3&         R,
                             const Vec3&         T,
                             const build::Model& o1,
                             const std::size_t   b1,
                             const build::Model& o2,
                             const std::size_t   b2,
                             const bool          rss_leaf_collision) const
{
    const build::BV& bv1 = o1.Bvs()[b1];
    const build::BV& bv2 = o2.Bvs()[b2];

    if (bv1.IsLeaf() && bv2.IsLeaf())
    {
        res.m_num_tri_tests++;

        if (rss_leaf_collision)
        {
            res.m_closer_than_tolerance = 1;
            res.m_distance = k_zero;
            return;
        }

        const Tri&                t1     = o1.OriginalTris()[bv1.TriIndex()];
        const Tri&                t2     = o2.OriginalTris()[bv2.TriIndex()];
        const std::array<Vec3, 3> t1_arr = {t1[0], t1[1], t1[2]};
        const std::array<Vec3, 3> t2_arr = {MxVpV(res.R(), t2[0], res.T()),
                                            MxVpV(res.R(), t2[1], res.T()),
                                            MxVpV(res.R(), t2[2], res.T())};
        Vec3                      p;
        Vec3                      q;
        Real                      tri_dist_sq = detail::TriDistSq(p, q, t1_arr, t2_arr);

        if (tri_dist_sq <= res.m_tolerance * res.m_tolerance)
        {
            res.m_closer_than_tolerance = 1;
            res.m_distance              = std::sqrt(tri_dist_sq);
            res.m_p1                    = p;
            res.m_p2                    = q;
        }
        return;
    }

    std::size_t a1;
    std::size_t a2;
    std::size_t c1;
    std::size_t c2;
    Mat3        R1;
    Mat3        R2;
    Vec3        T1;
    Vec3        T2;

    const Real sz1 = bv1.Size();
    const Real sz2 = bv2.Size();

    if (bv2.IsLeaf() || (!bv1.IsLeaf() && (sz1 > sz2)))
    {
        a2 = c2 = b2;
        a1      = bv1.FirstChild();
        c1      = bv1.FirstChild() + 1;

        const build::BV& ca1 = o1.Bvs()[a1];
        const build::BV& cc1 = o1.Bvs()[c1];

        R1 = MTxM(ca1.R(), R);
        T1 = MTxVmV(ca1.R(), T, ca1.Tr());

        R2 = MTxM(cc1.R(), R);
        T2 = MTxVmV(cc1.R(), T, cc1.Tr());
    }
    else
    {
        a1 = c1 = b1;
        a2      = bv2.FirstChild();
        c2      = bv2.FirstChild() + 1;

        const build::BV& ca2 = o2.Bvs()[a2];
        const build::BV& cc2 = o2.Bvs()[c2];

        R1 = MxM(R, ca2.R());
        T1 = MxVpV(R, ca2.Tr(), T);

        R2 = MxM(R, cc2.R());
        T2 = MxVpV(R, cc2.Tr(), T);
    }

    res.m_num_bv_tests += 2;

    const build::BV& bva1 = o1.Bvs()[a1];
    const build::BV& bva2 = o2.Bvs()[a2];
    const build::BV& bvc1 = o1.Bvs()[c1];
    const build::BV& bvc2 = o2.Bvs()[c2];

    const Real d1               = detail::RectDistSq(R1, T1, bva1.L(), bva2.L());
    const Real d2               = detail::RectDistSq(R2, T2, bvc1.L(), bvc2.L());
    Real       tolerance_plus_r = k_zero;

    if (d2 < d1)
    {
        tolerance_plus_r = res.m_tolerance + bvc1.Radius() + bvc2.Radius();
        if (d2 <= tolerance_plus_r * tolerance_plus_r)
        {
            QueryRecurse(res, R2, T2, o1, c1, o2, c2, rss_leaf_collision);
        }
        if (res.m_closer_than_tolerance)
        {
            return;
        }
        tolerance_plus_r = res.m_tolerance + bva1.Radius() + bva2.Radius();
        if (d1 <= tolerance_plus_r * tolerance_plus_r)
        {
            QueryRecurse(res, R1, T1, o1, a1, o2, a2, rss_leaf_collision);
        }
    }
    else
    {
        tolerance_plus_r = res.m_tolerance + bva1.Radius() + bva2.Radius();
        if (d1 <= tolerance_plus_r * tolerance_plus_r)
        {
            QueryRecurse(res, R1, T1, o1, a1, o2, a2, rss_leaf_collision);
        }
        if (res.m_closer_than_tolerance)
        {
            return;
        }
        tolerance_plus_r = res.m_tolerance + bvc1.Radius() + bvc2.Radius();
        if (d2 <= tolerance_plus_r * tolerance_plus_r)
        {
            QueryRecurse(res, R2, T2, o1, c1, o2, c2, rss_leaf_collision);
        }
    }
}

} // namespace pqss
