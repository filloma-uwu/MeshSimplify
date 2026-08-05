#include <algorithm>
#include <array>
#include <cmath>
#include <limits>
#include <numeric>
#include <numbers>
#include <vector>

#include "pqss/model_pool.hpp"
#include "pqss_proxy_mesh/proxy_bvh_build.hpp"
#include "pqss/utils/mat_vec.hpp"

namespace pqss::build
{
constexpr Real k_real_max    = std::numeric_limits<Real>::max();
constexpr Real k_real_lowest = std::numeric_limits<Real>::lowest();

struct SahBucket
{
    int  count = 0;
    Vec3 minv  = {k_real_max, k_real_max, k_real_max};
    Vec3 maxv  = {k_real_lowest, k_real_lowest, k_real_lowest};

    void add(const Tri& t)
    {
        count++;
        minv[0] = std::min({minv[0], t.p1[0], t.p2[0], t.p3[0]});
        maxv[0] = std::max({maxv[0], t.p1[0], t.p2[0], t.p3[0]});
        minv[1] = std::min({minv[1], t.p1[1], t.p2[1], t.p3[1]});
        maxv[1] = std::max({maxv[1], t.p1[1], t.p2[1], t.p3[1]});
        minv[2] = std::min({minv[2], t.p1[2], t.p2[2], t.p3[2]});
        maxv[2] = std::max({maxv[2], t.p1[2], t.p2[2], t.p3[2]});
    }
};

static Real ComputeArea(const Vec3& minv, const Vec3& maxv)
{
    Vec3 d = maxv - minv;
    if (d[0] < 0 || d[1] < 0 || d[2] < 0)
    {
        return 0;
    }
    return static_cast<Real>(2.0 * (d[0] * d[1] + d[1] * d[2] + d[2] * d[0]) + 1e-9);
}

static void MergeBounds(Vec3& res_min, Vec3& res_max, const Vec3& minv, const Vec3& maxv)
{
    res_min[0] = std::min(res_min[0], minv[0]);
    res_min[1] = std::min(res_min[1], minv[1]);
    res_min[2] = std::min(res_min[2], minv[2]);

    res_max[0] = std::max(res_max[0], maxv[0]);
    res_max[1] = std::max(res_max[1], maxv[1]);
    res_max[2] = std::max(res_max[2], maxv[2]);
}

Mat3 Model::GetCovarianceTriverts(const std::vector<SubTri>& tri_source,
                                  const std::size_t          first_tri,
                                  const std::size_t          num_tris) const
{
    Vec3 S1 = {0, 0, 0};
    Mat3 S2 = {{{{0, 0, 0}}, {{0, 0, 0}}, {{0, 0, 0}}}};

    for (int i = 0; i < num_tris; ++i)
    {
        const SubTri&             t = tri_source[first_tri + i];
        const std::array<Vec3, 3> p = {t.p1, t.p2, t.p3};
        for (int j = 0; j < 3; ++j)
        {
            S1[j] += p[0][j] + p[1][j] + p[2][j];
            for (int k = 0; k < 3; ++k)
            {
                S2[j][k] += (p[0][j] * p[0][k] + p[1][j] * p[1][k] + p[2][j] * p[2][k]);
            }
        }
    }

    Mat3 M{};
    Real n = static_cast<Real>(3 * num_tris);
    for (int i = 0; i < 3; ++i)
    {
        for (int j = 0; j < 3; ++j)
        {
            M[i][j] = S2[i][j] - (S1[i] * S1[j]) / n;
        }
    }
    return M;
}

Mat3 Model::GetCovarianceTriverts(const std::size_t first_tri, const std::size_t num_tris) const
{
    return GetCovarianceTriverts(m_tris, first_tri, num_tris);
}

void Model::BuildRecurse(const std::size_t bn, const std::size_t first_tri, const std::size_t num_tris)
{
    Mat3 C      = GetCovarianceTriverts(first_tri, num_tris);
    auto [E, s] = Meigen(C);

    int min_idx = 0, mid_idx = 1, max_idx = 2;
    if (s[max_idx] < s[mid_idx])
        std::swap(max_idx, mid_idx);
    if (s[max_idx] < s[min_idx])
        std::swap(max_idx, min_idx);
    if (s[mid_idx] < s[min_idx])
        std::swap(mid_idx, min_idx);

    Mat3 R{};
    for (int i = 0; i < 3; ++i)
    {
        R[i][0] = E[i][max_idx];
        R[i][1] = E[i][mid_idx];
    }
    R[0][2] = R[1][0] * R[2][1] - R[2][0] * R[1][1];
    R[1][2] = R[2][0] * R[0][1] - R[0][0] * R[2][1];
    R[2][2] = R[0][0] * R[1][1] - R[1][0] * R[0][1];

    this->m_b[bn].FitToTris(R, this->m_tris, first_tri, num_tris, m_build_strategy);

    if (num_tris <= 1)
    {
        this->m_b[bn].SetLeaf(this->m_tris[first_tri].original_tri_id);
        if (bn >= this->m_leaf_subtri_indices.size())
        {
            this->m_leaf_subtri_indices.resize(bn + 1, std::numeric_limits<std::size_t>::max());
        }
        this->m_leaf_subtri_indices[bn] = first_tri;
        return;
    }

    constexpr int k_num_bins = 16;
    Real          best_cost  = std::numeric_limits<Real>::max();
    Vec3          best_axis{};
    int           best_split_bin = -1;
    Real          best_min_c = 0, best_max_c = 0;

    const std::array<Vec3, 6> split_axes = {
        Vec3{1, 0, 0},
        Vec3{0, 1, 0},
        Vec3{0, 0, 1},
        Vec3{R[0][0], R[1][0], R[2][0]},
        Vec3{R[0][1], R[1][1], R[2][1]},
        Vec3{R[0][2], R[1][2], R[2][2]},
    };
    const bool proxy_build = pqss_proxy_mesh::proxyBvhBuildEnabled();
    const int axis_count = proxy_build ? static_cast<int>(split_axes.size()) : 3;

    for (int axis_index = 0; axis_index < axis_count; ++axis_index)
    {
        const Vec3& axis = split_axes[axis_index];
        Real min_c = k_real_max;
        Real max_c = k_real_lowest;

        for (int i = 0; i < num_tris; ++i)
        {
            const SubTri& tri = this->m_tris[first_tri + i];
            const Vec3 centroid = (tri.p1 + tri.p2 + tri.p3) / Real(3);
            const Real cent = VdotV(centroid, axis);
            min_c     = std::min(min_c, cent);
            max_c     = std::max(max_c, cent);
        }

        if (max_c - min_c < 1e-7)
            continue;

        std::array<SahBucket, k_num_bins> bins;
        for (int i = 0; i < num_tris; ++i)
        {
            const SubTri& tri = this->m_tris[first_tri + i];
            const Vec3 centroid = (tri.p1 + tri.p2 + tri.p3) / Real(3);
            const Real cent = VdotV(centroid, axis);
            int  b_idx = static_cast<int>((k_num_bins - 1) * (cent - min_c) / (max_c - min_c));
            b_idx      = std::clamp(b_idx, 0, k_num_bins - 1);
            bins[b_idx].add(this->m_tris[first_tri + i]);
        }

        std::array<Vec3, k_num_bins> pre_min{}, pre_max{};
        std::array<int, k_num_bins>  pre_cnt{};
        Vec3                         cur_min = {k_real_max, k_real_max, k_real_max};
        Vec3                         cur_max = {k_real_lowest, k_real_lowest, k_real_lowest};
        int                          cur_cnt = 0;
        for (int i = 0; i < k_num_bins; ++i)
        {
            if (bins[i].count > 0)
                MergeBounds(cur_min, cur_max, bins[i].minv, bins[i].maxv);
            cur_cnt += bins[i].count;
            pre_min[i] = cur_min;
            pre_max[i] = cur_max;
            pre_cnt[i] = cur_cnt;
        }

        int  p_cnt    = 0;
        Vec3 post_min = {k_real_max, k_real_max, k_real_max};
        Vec3 post_max = {k_real_lowest, k_real_lowest, k_real_lowest};
        for (std::size_t i = k_num_bins - 1; i > 0; --i)
        {
            if (bins[i].count > 0)
                MergeBounds(post_min, post_max, bins[i].minv, bins[i].maxv);
            p_cnt += bins[i].count;
            if (pre_cnt[i - 1] == 0 || p_cnt == 0)
                continue;
            Real cost =
                pre_cnt[i - 1] * ComputeArea(pre_min[i - 1], pre_max[i - 1]) + p_cnt * ComputeArea(post_min, post_max);
            if (cost < best_cost)
            {
                best_cost      = cost;
                best_axis      = axis;
                best_split_bin = static_cast<int>(i);
                best_min_c     = min_c;
                best_max_c     = max_c;
            }
        }
    }

    std::size_t split_idx;
    auto        start_it = this->m_tris.begin() + first_tri;
    auto        end_it   = start_it + num_tris;

    if (best_split_bin == -1)
    {
        split_idx = num_tris / 2;
    }
    else
    {
        auto it = std::partition(
            start_it,
            end_it,
            [&](const Tri& t)
            {
                const Vec3 centroid = (t.p1 + t.p2 + t.p3) / Real(3);
                const Real cent = VdotV(centroid, best_axis);
                int  b_idx = static_cast<int>((k_num_bins - 1) * (cent - best_min_c) / (best_max_c - best_min_c));
                return b_idx < best_split_bin;
            });
        split_idx = static_cast<int>(std::distance(start_it, it));
    }

    if (split_idx <= 0 || split_idx >= num_tris)
        split_idx = num_tris / 2;

    std::size_t lc = m_next_bv_idx;
    m_next_bv_idx += 2;

    this->m_b[bn].SetInternal(lc);

    BuildRecurse(lc, first_tri, split_idx);
    BuildRecurse(lc + 1, first_tri + split_idx, num_tris - split_idx);
}

void Model::FitB12vNode(const std::vector<SubTri>& ordered_tris,
                        const std::size_t          first_tri,
                        const std::size_t          num_tris,
                        const std::size_t          node)
{
    Mat3 R{};
    if (num_tris <= 2)
    {
        R = Meident();
    }
    else
    {
        Mat3 C      = GetCovarianceTriverts(ordered_tris, first_tri, num_tris);
        auto [E, s] = Meigen(C);

        int min_idx = 0, mid_idx = 1, max_idx = 2;
        if (s[max_idx] < s[mid_idx])
            std::swap(max_idx, mid_idx);
        if (s[max_idx] < s[min_idx])
            std::swap(max_idx, min_idx);
        if (s[mid_idx] < s[min_idx])
            std::swap(mid_idx, min_idx);

        for (int i = 0; i < 3; ++i)
        {
            R[i][0] = E[i][max_idx];
            R[i][1] = E[i][mid_idx];
        }
        R[0][2] = R[1][0] * R[2][1] - R[2][0] * R[1][1];
        R[1][2] = R[2][0] * R[0][1] - R[0][0] * R[2][1];
        R[2][2] = R[0][0] * R[1][1] - R[1][0] * R[0][1];
    }

    BV& bv = m_b12v_rss_b[node];
    bv.FitToTris(R, ordered_tris, first_tri, num_tris, BuildStrategy::Fast);
    bv.m_first_child = -1;
    bv.m_tri0        = -1;
    bv.m_tri1        = -1;
}

std::pair<std::size_t, std::size_t> Model::BuildB12vNodeRanges(const std::vector<SubTri>&   ordered_tris,
                                                               const std::vector<B12vLeaf>& leaves,
                                                               const std::size_t            node,
                                                               const std::size_t            leaf_count)
{
    const std::size_t first_leaf = leaf_count - 1;
    if (node >= first_leaf)
    {
        const std::size_t leaf_index = node - first_leaf;
        const B12vLeaf&   leaf       = leaves[leaf_index];
        const std::size_t first_tri  = leaf.first_tri;
        const std::size_t count      = leaf.count;
        FitB12vNode(ordered_tris, first_tri, count, node);
        BV& bv    = m_b12v_rss_b[node];
        bv.m_tri0 = static_cast<int>(ordered_tris[first_tri].original_tri_id);
        bv.m_tri1 = count > 1 ? static_cast<int>(ordered_tris[first_tri + 1].original_tri_id) : -1;
        return {first_tri, count};
    }

    const auto [left_first, left_count]   = BuildB12vNodeRanges(ordered_tris, leaves, node * 2 + 1, leaf_count);
    const auto [right_first, right_count] = BuildB12vNodeRanges(ordered_tris, leaves, node * 2 + 2, leaf_count);
    const std::size_t first_tri           = std::min(left_first, right_first);
    const std::size_t end_tri             = std::max(left_first + left_count, right_first + right_count);
    FitB12vNode(ordered_tris, first_tri, end_tri - first_tri, node);
    return {first_tri, end_tri - first_tri};
}

Return Model::BuildB12vRss()
{
    m_b12v_rss_b.clear();
    m_b12v_rss_first_single_leaf = 0;

    if (m_tris.empty())
        return Return::EmptyModel;

    std::size_t leaf_count = 1;
    while ((leaf_count << 1) <= m_tris.size())
        leaf_count <<= 1;

    const std::size_t node_count = leaf_count * 2 - 1;
    try
    {
        m_b12v_rss_b.resize(node_count);
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    std::vector<SubTri> ordered_tris = m_tris;
    Vec3                min_c        = {k_real_max, k_real_max, k_real_max};
    Vec3                max_c        = {k_real_lowest, k_real_lowest, k_real_lowest};
    for (const SubTri& tri : ordered_tris)
    {
        const Vec3 c = {(tri.p1[0] + tri.p2[0] + tri.p3[0]) / Real(3),
                        (tri.p1[1] + tri.p2[1] + tri.p3[1]) / Real(3),
                        (tri.p1[2] + tri.p2[2] + tri.p3[2]) / Real(3)};
        min_c[0]     = std::min(min_c[0], c[0]);
        min_c[1]     = std::min(min_c[1], c[1]);
        min_c[2]     = std::min(min_c[2], c[2]);
        max_c[0]     = std::max(max_c[0], c[0]);
        max_c[1]     = std::max(max_c[1], c[1]);
        max_c[2]     = std::max(max_c[2], c[2]);
    }
    auto expand_bits = [](std::uint32_t v)
    {
        v = (v | (v << 16)) & 0x030000FFu;
        v = (v | (v << 8)) & 0x0300F00Fu;
        v = (v | (v << 4)) & 0x030C30C3u;
        v = (v | (v << 2)) & 0x09249249u;
        return v;
    };
    auto morton_code = [&](const SubTri& tri)
    {
        const Vec3 c        = {(tri.p1[0] + tri.p2[0] + tri.p3[0]) / Real(3),
                               (tri.p1[1] + tri.p2[1] + tri.p3[1]) / Real(3),
                               (tri.p1[2] + tri.p2[2] + tri.p3[2]) / Real(3)};
        auto       quantize = [](const Real value, const Real min_value, const Real max_value)
        {
            if (max_value <= min_value)
                return std::uint32_t{0};
            const Real normalized = std::clamp((value - min_value) / (max_value - min_value), Real(0), Real(1));
            return static_cast<std::uint32_t>(normalized * Real(1023));
        };
        const std::uint32_t x = quantize(c[0], min_c[0], max_c[0]);
        const std::uint32_t y = quantize(c[1], min_c[1], max_c[1]);
        const std::uint32_t z = quantize(c[2], min_c[2], max_c[2]);
        return (expand_bits(x) << 2) | (expand_bits(y) << 1) | expand_bits(z);
    };
    std::stable_sort(ordered_tris.begin(),
                     ordered_tris.end(),
                     [&](const SubTri& lhs, const SubTri& rhs)
                     {
                         const std::uint32_t lm = morton_code(lhs);
                         const std::uint32_t rm = morton_code(rhs);
                         if (lm != rm)
                             return lm < rm;
                         return lhs.original_tri_id < rhs.original_tri_id;
                     });

    std::vector<B12vLeaf> leaves;
    try
    {
        leaves.reserve(leaf_count);
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    try
    {
        const std::size_t merge_count = ordered_tris.size() - leaf_count;
        std::size_t       tri_i       = 0;
        for (std::size_t leaf_i = 0; leaf_i < leaf_count; ++leaf_i)
        {
            const std::size_t prev_merges = leaf_i * merge_count / leaf_count;
            const std::size_t next_merges = (leaf_i + 1) * merge_count / leaf_count;
            const std::size_t count       = next_merges > prev_merges ? 2 : 1;
            leaves.push_back({tri_i, count});
            tri_i += count;
        }
        if (tri_i != ordered_tris.size())
            return Return::UnknownError;
    }
    catch (const std::bad_alloc&)
    {
        return Return::OutOfMemory;
    }

    if (leaves.size() != leaf_count)
        return Return::UnknownError;

    BuildB12vNodeRanges(ordered_tris, leaves, 0, leaf_count);
    m_b12v_rss_first_single_leaf = m_b12v_rss_b.size();
    return Return::Success;
}

void Model::MakeParentRelative(const std::size_t bn, const Mat3& parentR, const Vec3& parentTr)
{
    BV& bv = this->m_b[bn];
    if (!bv.IsLeaf())
    {
        MakeParentRelative(bv.FirstChild(), bv.m_R, bv.m_tr);
        MakeParentRelative(bv.FirstChild() + 1, bv.m_R, bv.m_tr);
    }

    Mat3 R_rel = MTxM(parentR, bv.m_R);
    Vec3 T_rel = MTxV(parentR, bv.m_tr - parentTr);

    bv.m_R  = R_rel;
    bv.m_tr = T_rel;
}

Return Model::BuildBV()
{
    if (m_tris.empty())
        return Return::EmptyModel;

    m_b.clear();
    m_b.resize(2 * m_tris.size() - 1);
    m_leaf_subtri_indices.clear();
    m_leaf_subtri_indices.resize(m_b.size(), std::numeric_limits<std::size_t>::max());

    m_next_bv_idx = 1;

    BuildRecurse(0, 0, static_cast<std::size_t>(m_tris.size()));

    MakeParentRelative(0, Meident(), Veident());

    return Return::Success;
}

void BV::FitToTrisFast(const Mat3&                O,
                       const std::vector<SubTri>& all_tris,
                       const std::size_t          first_tri,
                       const std::size_t          num_tris)
{
    this->m_R = O;

    const std::size_t num_points = 3 * num_tris;
    std::vector<Vec3> P(num_points);

    for (std::size_t i = 0; i < num_tris; ++i)
    {
        const Tri& t = all_tris[first_tri + i];
        P[i * 3 + 0] = MTxV(m_R, t.p1);
        P[i * 3 + 1] = MTxV(m_R, t.p2);
        P[i * 3 + 2] = MTxV(m_R, t.p3);
    }

    Real minx, maxx, miny, maxy, minz, maxz;
    Vec3 c{};

    Real cz, radsqr;
    minz = maxz = P[0][2];
    for (std::size_t i = 1; i < num_points; ++i)
    {
        if (P[i][2] < minz)
            minz = P[i][2];
        else if (P[i][2] > maxz)
            maxz = P[i][2];
    }

    this->m_r = static_cast<Real>(0.5) * (maxz - minz);
    radsqr    = m_r * m_r;
    cz        = static_cast<Real>(0.5) * (maxz + minz);

    std::size_t minindex = 0, maxindex = 0;
    for (std::size_t i = 1; i < num_points; ++i)
    {
        if (P[i][0] < P[minindex][0])
            minindex = i;
        else if (P[i][0] > P[maxindex][0])
            maxindex = i;
    }

    Real x, dz;
    dz   = P[minindex][2] - cz;
    minx = P[minindex][0] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
    dz   = P[maxindex][2] - cz;
    maxx = P[maxindex][0] - std::sqrt(std::max(radsqr - dz * dz, k_zero));

    for (std::size_t i = 0; i < num_points; ++i)
    {
        if (P[i][0] < minx)
        {
            dz = P[i][2] - cz;
            x  = P[i][0] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
            if (x < minx)
                minx = x;
        }
        if (P[i][0] > maxx)
        {
            dz = P[i][2] - cz;
            x  = P[i][0] - std::sqrt(std::max(radsqr - dz * dz, k_zero));
            if (x > maxx)
                maxx = x;
        }
    }

    minindex = maxindex = 0;
    for (std::size_t i = 1; i < num_points; ++i)
    {
        if (P[i][1] < P[minindex][1])
            minindex = i;
        else if (P[i][1] > P[maxindex][1])
            maxindex = i;
    }

    Real y;
    dz   = P[minindex][2] - cz;
    miny = P[minindex][1] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
    dz   = P[maxindex][2] - cz;
    maxy = P[maxindex][1] - std::sqrt(std::max(radsqr - dz * dz, k_zero));

    for (std::size_t i = 0; i < num_points; ++i)
    {
        if (P[i][1] < miny)
        {
            dz = P[i][2] - cz;
            y  = P[i][1] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
            if (y < miny)
                miny = y;
        }
        if (P[i][1] > maxy)
        {
            dz = P[i][2] - cz;
            y  = P[i][1] - std::sqrt(std::max(radsqr - dz * dz, k_zero));
            if (y > maxy)
                maxy = y;
        }
    }

    Real       dx, dy, u, t;
    const Real a_val = std::sqrt(static_cast<Real>(0.5));
    for (std::size_t i = 0; i < num_points; ++i)
    {
        if (P[i][0] > maxx)
        {
            if (P[i][1] > maxy)
            {
                dx = P[i][0] - maxx;
                dy = P[i][1] - maxy;
                u  = dx * a_val + dy * a_val;
                t  = (a_val * u - dx) * (a_val * u - dx) + (a_val * u - dy) * (a_val * u - dy) +
                     (cz - P[i][2]) * (cz - P[i][2]);
                u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                if (u > 0)
                {
                    maxx += u * a_val;
                    maxy += u * a_val;
                }
            }
            else if (P[i][1] < miny)
            {
                dx = P[i][0] - maxx;
                dy = P[i][1] - miny;
                u  = dx * a_val - dy * a_val;
                t  = (a_val * u - dx) * (a_val * u - dx) + (-a_val * u - dy) * (-a_val * u - dy) +
                     (cz - P[i][2]) * (cz - P[i][2]);
                u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                if (u > 0)
                {
                    maxx += u * a_val;
                    miny -= u * a_val;
                }
            }
        }
        else if (P[i][0] < minx)
        {
            if (P[i][1] > maxy)
            {
                dx = P[i][0] - minx;
                dy = P[i][1] - maxy;
                u  = dy * a_val - dx * a_val;
                t  = (-a_val * u - dx) * (-a_val * u - dx) + (a_val * u - dy) * (a_val * u - dy) +
                     (cz - P[i][2]) * (cz - P[i][2]);
                u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                if (u > 0)
                {
                    minx -= u * a_val;
                    maxy += u * a_val;
                }
            }
            else if (P[i][1] < miny)
            {
                dx = P[i][0] - minx;
                dy = P[i][1] - miny;
                u  = -dx * a_val - dy * a_val;
                t  = (-a_val * u - dx) * (-a_val * u - dx) + (-a_val * u - dy) * (-a_val * u - dy) +
                     (cz - P[i][2]) * (cz - P[i][2]);
                u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                if (u > 0)
                {
                    minx -= u * a_val;
                    miny -= u * a_val;
                }
            }
        }
    }

    c[0]       = minx;
    c[1]       = miny;
    c[2]       = cz;
    this->m_tr = MxV(m_R, c);

    this->m_l[0] = std::max(maxx - minx, k_zero);
    this->m_l[1] = std::max(maxy - miny, k_zero);

    this->m_size = CalcSize();
}

void BV::FitToTrisOptimized(const Mat3&                O,
                            const std::vector<SubTri>& all_tris,
                            const std::size_t          first_tri,
                            const std::size_t          num_tris)
{
    const std::size_t num_points = 3 * num_tris;
    std::vector<Vec3> points(num_points);

    for (std::size_t i = 0; i < num_tris; ++i)
    {
        const Tri& t      = all_tris[first_tri + i];
        points[i * 3 + 0] = t.p1;
        points[i * 3 + 1] = t.p2;
        points[i * 3 + 2] = t.p3;
    }

    auto rotate_local = [](const Mat3& frame, const int axis, const Real angle)
    {
        const Real cs  = std::cos(angle);
        const Real sn  = std::sin(angle);
        Mat3       out = frame;
        const int  a   = (axis + 1) % 3;
        const int  b   = (axis + 2) % 3;
        for (int row = 0; row < 3; ++row)
        {
            out[row][a] = frame[row][a] * cs + frame[row][b] * sn;
            out[row][b] = frame[row][b] * cs - frame[row][a] * sn;
        }
        return out;
    };

    struct Candidate
    {
        Mat3                R{};
        Vec3                Tr{};
        std::array<Real, 2> l{};
        Real                r    = k_zero;
        Real                size = std::numeric_limits<Real>::max();
    };

    auto evaluate_frame = [&](const Mat3& frame, Candidate& best_candidate)
    {
        std::vector<Vec3> P(num_points);
        for (std::size_t i = 0; i < num_points; ++i)
        {
            P[i] = MTxV(frame, points[i]);
        }

        Real minz, maxz;
        minz = maxz = P[0][2];
        for (std::size_t i = 1; i < num_points; ++i)
        {
            if (P[i][2] < minz)
                minz = P[i][2];
            else if (P[i][2] > maxz)
                maxz = P[i][2];
        }

        const Real min_r = static_cast<Real>(0.5) * (maxz - minz);
        const Real mid_z = static_cast<Real>(0.5) * (maxz + minz);

        auto evaluate_radius = [&](const Real r, const Real cz, Candidate& radius_best)
        {
            Real       minx, maxx, miny, maxy;
            Vec3       c{};
            const Real radsqr = r * r;

            std::size_t minindex = 0, maxindex = 0;
            for (std::size_t i = 1; i < num_points; ++i)
            {
                if (P[i][0] < P[minindex][0])
                    minindex = i;
                else if (P[i][0] > P[maxindex][0])
                    maxindex = i;
            }

            Real x, dz;
            dz   = P[minindex][2] - cz;
            minx = P[minindex][0] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
            dz   = P[maxindex][2] - cz;
            maxx = P[maxindex][0] - std::sqrt(std::max(radsqr - dz * dz, k_zero));

            for (std::size_t i = 0; i < num_points; ++i)
            {
                if (P[i][0] < minx)
                {
                    dz = P[i][2] - cz;
                    x  = P[i][0] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
                    if (x < minx)
                        minx = x;
                }
                if (P[i][0] > maxx)
                {
                    dz = P[i][2] - cz;
                    x  = P[i][0] - std::sqrt(std::max(radsqr - dz * dz, k_zero));
                    if (x > maxx)
                        maxx = x;
                }
            }

            minindex = maxindex = 0;
            for (std::size_t i = 1; i < num_points; ++i)
            {
                if (P[i][1] < P[minindex][1])
                    minindex = i;
                else if (P[i][1] > P[maxindex][1])
                    maxindex = i;
            }

            Real y;
            dz   = P[minindex][2] - cz;
            miny = P[minindex][1] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
            dz   = P[maxindex][2] - cz;
            maxy = P[maxindex][1] - std::sqrt(std::max(radsqr - dz * dz, k_zero));

            for (std::size_t i = 0; i < num_points; ++i)
            {
                if (P[i][1] < miny)
                {
                    dz = P[i][2] - cz;
                    y  = P[i][1] + std::sqrt(std::max(radsqr - dz * dz, k_zero));
                    if (y < miny)
                        miny = y;
                }
                if (P[i][1] > maxy)
                {
                    dz = P[i][2] - cz;
                    y  = P[i][1] - std::sqrt(std::max(radsqr - dz * dz, k_zero));
                    if (y > maxy)
                        maxy = y;
                }
            }

            // Overlapping side constraints mean that the core segment should
            // collapse to a point before points are classified by corner.
            auto collapse_overlapping_interval = [](Real& minimum, Real& maximum)
            {
                if (minimum > maximum)
                {
                    const Real midpoint = std::midpoint(minimum, maximum);
                    minimum             = midpoint;
                    maximum             = midpoint;
                }
            };
            collapse_overlapping_interval(minx, maxx);
            collapse_overlapping_interval(miny, maxy);

            Real       dx, dy, u, t;
            const Real a_val = std::sqrt(static_cast<Real>(0.5));
            for (std::size_t i = 0; i < num_points; ++i)
            {
                if (P[i][0] > maxx)
                {
                    if (P[i][1] > maxy)
                    {
                        dx = P[i][0] - maxx;
                        dy = P[i][1] - maxy;
                        u  = dx * a_val + dy * a_val;
                        t  = (a_val * u - dx) * (a_val * u - dx) + (a_val * u - dy) * (a_val * u - dy) +
                             (cz - P[i][2]) * (cz - P[i][2]);
                        u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                        if (u > 0)
                        {
                            maxx += u * a_val;
                            maxy += u * a_val;
                        }
                    }
                    else if (P[i][1] < miny)
                    {
                        dx = P[i][0] - maxx;
                        dy = P[i][1] - miny;
                        u  = dx * a_val - dy * a_val;
                        t  = (a_val * u - dx) * (a_val * u - dx) + (-a_val * u - dy) * (-a_val * u - dy) +
                             (cz - P[i][2]) * (cz - P[i][2]);
                        u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                        if (u > 0)
                        {
                            maxx += u * a_val;
                            miny -= u * a_val;
                        }
                    }
                }
                else if (P[i][0] < minx)
                {
                    if (P[i][1] > maxy)
                    {
                        dx = P[i][0] - minx;
                        dy = P[i][1] - maxy;
                        u  = dy * a_val - dx * a_val;
                        t  = (-a_val * u - dx) * (-a_val * u - dx) + (a_val * u - dy) * (a_val * u - dy) +
                             (cz - P[i][2]) * (cz - P[i][2]);
                        u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                        if (u > 0)
                        {
                            minx -= u * a_val;
                            maxy += u * a_val;
                        }
                    }
                    else if (P[i][1] < miny)
                    {
                        dx = P[i][0] - minx;
                        dy = P[i][1] - miny;
                        u  = -dx * a_val - dy * a_val;
                        t  = (-a_val * u - dx) * (-a_val * u - dx) + (-a_val * u - dy) * (-a_val * u - dy) +
                             (cz - P[i][2]) * (cz - P[i][2]);
                        u  = u - std::sqrt(std::max(radsqr - t, k_zero));
                        if (u > 0)
                        {
                            minx -= u * a_val;
                            miny -= u * a_val;
                        }
                    }
                }
            }

            Candidate candidate;
            candidate.R    = frame;
            c[0]           = minx;
            c[1]           = miny;
            c[2]           = cz;
            candidate.Tr   = MxV(candidate.R, c);
            candidate.l[0] = std::max(maxx - minx, k_zero);
            candidate.l[1] = std::max(maxy - miny, k_zero);
            candidate.r    = r;
            candidate.size = candidate.l[0] * candidate.l[1] +
                             std::numbers::pi_v<Real> * candidate.r *
                                 (candidate.l[0] + candidate.l[1] + static_cast<Real>(2) * candidate.r);

            if (candidate.size < radius_best.size)
            {
                radius_best = candidate;
            }
        };

        Real max_plane_dist_sq = k_zero;
        for (const Vec3& p : P)
        {
            const Real dx     = p[0];
            const Real dy     = p[1];
            const Real dz     = p[2] - mid_z;
            max_plane_dist_sq = std::max(max_plane_dist_sq, dx * dx + dy * dy + dz * dz);
        }

        const Real max_r = std::sqrt(max_plane_dist_sq);
        Candidate  radius_best;
        evaluate_radius(min_r, mid_z, radius_best);
        evaluate_radius(max_r, mid_z, radius_best);

        if (max_r > min_r + static_cast<Real>(1e-12))
        {
            Real           lo      = min_r;
            Real           hi      = max_r;
            constexpr Real inv_phi = static_cast<Real>(0.6180339887498948482);
            Real           c       = hi - (hi - lo) * inv_phi;
            Real           d       = lo + (hi - lo) * inv_phi;
            Candidate      c_candidate;
            Candidate      d_candidate;
            evaluate_radius(c, mid_z, c_candidate);
            evaluate_radius(d, mid_z, d_candidate);

            for (int iter = 0; iter < 16; ++iter)
            {
                if (c_candidate.size < d_candidate.size)
                {
                    hi          = d;
                    d           = c;
                    d_candidate = c_candidate;
                    c           = hi - (hi - lo) * inv_phi;
                    c_candidate = Candidate{};
                    evaluate_radius(c, mid_z, c_candidate);
                }
                else
                {
                    lo          = c;
                    c           = d;
                    c_candidate = d_candidate;
                    d           = lo + (hi - lo) * inv_phi;
                    d_candidate = Candidate{};
                    evaluate_radius(d, mid_z, d_candidate);
                }
            }

            if (c_candidate.size < radius_best.size)
            {
                radius_best = c_candidate;
            }
            if (d_candidate.size < radius_best.size)
            {
                radius_best = d_candidate;
            }
        }

        if (radius_best.size < best_candidate.size)
        {
            best_candidate = radius_best;
        }
    };

    Candidate best;
    evaluate_frame(O, best);

    auto optimize_rotation_line = [&](const int axis, const Real radius)
    {
        Real           lo      = -radius;
        Real           hi      = radius;
        constexpr Real inv_phi = static_cast<Real>(0.6180339887498948482);
        Real           c       = hi - (hi - lo) * inv_phi;
        Real           d       = lo + (hi - lo) * inv_phi;

        Candidate c_candidate;
        Candidate d_candidate;
        evaluate_frame(rotate_local(best.R, axis, c), c_candidate);
        evaluate_frame(rotate_local(best.R, axis, d), d_candidate);

        for (int iter = 0; iter < 10; ++iter)
        {
            if (c_candidate.size < d_candidate.size)
            {
                hi          = d;
                d           = c;
                d_candidate = c_candidate;
                c           = hi - (hi - lo) * inv_phi;
                c_candidate = Candidate{};
                evaluate_frame(rotate_local(best.R, axis, c), c_candidate);
            }
            else
            {
                lo          = c;
                c           = d;
                c_candidate = d_candidate;
                d           = lo + (hi - lo) * inv_phi;
                d_candidate = Candidate{};
                evaluate_frame(rotate_local(best.R, axis, d), d_candidate);
            }
        }

        Candidate line_best = best;
        if (c_candidate.size < line_best.size)
        {
            line_best = c_candidate;
        }
        if (d_candidate.size < line_best.size)
        {
            line_best = d_candidate;
        }
        if (line_best.size < best.size)
        {
            best = line_best;
            return true;
        }
        return false;
    };

    const std::array<Real, 3> trust_radii = {static_cast<Real>(0.18),
                                             static_cast<Real>(0.09),
                                             static_cast<Real>(0.045)};
    for (const Real radius : trust_radii)
    {
        for (int pass = 0; pass < 2; ++pass)
        {
            bool improved = false;
            for (int axis = 0; axis < 3; ++axis)
            {
                improved = optimize_rotation_line(axis, radius) || improved;
            }
            if (!improved)
            {
                break;
            }
        }
    }

    this->m_R    = best.R;
    this->m_tr   = best.Tr;
    this->m_l    = best.l;
    this->m_r    = best.r;
    this->m_size = best.size;
}

void BV::FitToTris(const Mat3&                O,
                   const std::vector<SubTri>& all_tris,
                   const std::size_t          first_tri,
                   const std::size_t          num_tris,
                   const BuildStrategy        build_strategy)
{
    if (build_strategy == BuildStrategy::Fast)
    {
        FitToTrisFast(O, all_tris, first_tri, num_tris);
    }
    else
    {
        FitToTrisOptimized(O, all_tris, first_tri, num_tris);
    }
}

} // namespace pqss::build
