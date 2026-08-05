#pragma once

#include <algorithm>
#include <array>
#include <cassert>
#include <cstddef>
#include <numbers>
#include <vector>

#include "pqss/config/constant.hpp"

namespace pqss
{
class ModelPool;
}

namespace pqss::build
{

struct Tri
{
    Vec3 p1{};
    Vec3 p2{};
    Vec3 p3{};

    Tri() = default;
    Tri(const Vec3& first, const Vec3& second, const Vec3& third)
        : p1{first}, p2{second}, p3{third}
    {
    }

    const Vec3& operator[](const std::size_t index) const
    {
        assert(index < 3);
        if (index == 0) return p1;
        if (index == 1) return p2;
        return p3;
    }
};

struct SubTri : public Tri
{
    std::size_t original_tri_id = 0;

    SubTri() = default;
    SubTri(const Vec3& first,
           const Vec3& second,
           const Vec3& third,
           const std::size_t original_id)
        : Tri{first, second, third}, original_tri_id{original_id}
    {
    }
};

class BV
{
public:
    const Mat3& R() const
    {
        return m_R;
    }

    const Vec3& Tr() const
    {
        return m_tr;
    }

    const std::array<Real, 2>& L() const
    {
        return m_l;
    }

    Real Radius() const
    {
        return m_r;
    }

    Real Size() const
    {
        return m_size;
    }

    std::size_t FirstChild() const
    {
        return static_cast<std::size_t>(m_first_child);
    }

    std::size_t TriIndex() const
    {
        return static_cast<std::size_t>(-m_first_child - 1);
    }

    bool IsLeaf() const
    {
        return m_first_child < 0;
    }

    void SetInternal(const std::size_t first_child)
    {
        m_first_child = static_cast<int>(first_child);
    }

    void SetLeaf(const std::size_t triangle_index)
    {
        m_first_child = -static_cast<int>(triangle_index + 1);
    }

    int FirstChildRaw() const
    {
        return m_first_child;
    }

    int Tri0() const
    {
        return m_tri0;
    }

    int Tri1() const
    {
        return m_tri1;
    }

    bool Leaf() const
    {
        return IsLeaf();
    }

    // Project-local construction API. The implementation remains PQSS's
    // original Fast/Optimized RSS fitter in pqss_proxy_build.cpp.
    void FitToTris(const Mat3&                initial_frame,
                   const std::vector<SubTri>& triangles,
                   std::size_t                first_triangle,
                   std::size_t                triangle_count,
                   BuildStrategy              strategy);

    friend class Model;
    friend class ::pqss::ModelPool;

private:
    Mat3                m_R{};
    Vec3                m_tr{};
    std::array<Real, 2> m_l{};
    Real                m_r           = k_zero;
    Real                m_size        = k_zero;
    int                 m_first_child = 0;
    int                 m_tri0        = -1;
    int                 m_tri1        = -1;

    Real CalcSize() const
    {
        return m_l[0] * m_l[1] + std::numbers::pi_v<Real> * m_r * (m_l[0] + m_l[1] + 2 * m_r);
    }

    void FitToTrisFast(const Mat3&                initial_frame,
                       const std::vector<SubTri>& triangles,
                       std::size_t                first_triangle,
                       std::size_t                triangle_count);

    void FitToTrisOptimized(const Mat3&                initial_frame,
                            const std::vector<SubTri>& triangles,
                            std::size_t                first_triangle,
                            std::size_t                triangle_count);
};

} // namespace pqss::build
