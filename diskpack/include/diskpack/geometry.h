#pragma once

#include <boost/numeric/interval.hpp>
#include <memory>
#include <optional>
#include <vector>

namespace diskpack {

using BaseType = long double;

using Interval = boost::numeric::interval<
    BaseType,
    boost::numeric::interval_lib::policies<
        boost::numeric::interval_lib::save_state<
            boost::numeric::interval_lib::rounded_transc_exact<BaseType>>,
        boost::numeric::interval_lib::checking_base<BaseType>>>;

using IntervalPair = std::pair<Interval, Interval>;

const BaseType EPS = 1e-11;

const Interval ONE{1 - EPS, 1 + EPS};
const Interval ZERO{-EPS, EPS};


/// Spiral similarity operator
/// Used to compute coordinates of a disk tangent to two others

struct SpiralOp {

    Interval x;
    Interval y;

    SpiralOp(const Interval& x_, const Interval& y_);
    SpiralOp(const SpiralOp& other);
    SpiralOp(SpiralOp&& other);

    SpiralOp& operator=(const SpiralOp& other);
    SpiralOp& operator=(SpiralOp&& other);

    SpiralOp(
        const Interval& base_r,
        const Interval& prev_r,
        const Interval& next_r,
        const size_t& base_t,
        const size_t& prev_t,
        const size_t& next_t
    );

    IntervalPair operator*(const IntervalPair& vec) const noexcept;

    SpiralOp operator*(const SpiralOp& other) const noexcept;

    SpiralOp();
};



/// Disk representation

class Disk {

private:

    Interval center_x;
    Interval center_y;
    Interval radius;

    size_t disk_type;

public:

    Disk(
        Interval center_x_,
        Interval center_y_,
        Interval radius_,
        size_t disk_type_
    );

    Disk();

    Disk(const Disk& other) noexcept;
    Disk(Disk&& other) noexcept;

    Disk& operator=(const Disk& other);
    Disk& operator=(Disk&& other);

    Interval getNorm() const;

    const Interval& getRadius() const;
    const Interval& getCenterX() const;
    const Interval& getCenterY() const;

    size_t getType() const;

    bool intersects(const Disk& other) const;
    bool tangent(const Disk& other) const;
    bool disjoint(const Disk& other) const;

    BaseType precision() const;
};



using DiskPtr = std::shared_ptr<Disk>;

bool compareDiskNorm(const DiskPtr a, const DiskPtr b);


/// Region of radii search

class RadiusRegion {

private:

    std::vector<Interval> intervals;

public:

    RadiusRegion(const std::vector<Interval>& intervals_);
    RadiusRegion(std::vector<Interval>&& intervals_);

    // ===== НОВЫЙ API (твой текущий) =====
    const std::vector<Interval>& getIntervals() const;

    bool isNarrowEnough(BaseType lower_bound) const;
    bool isTooWide(BaseType upper_bound) const;

    Interval getMinInterval() const;
    Interval getMaxInterval() const;

    void split(
        std::vector<RadiusRegion>& regions,
        size_t k,
        std::optional<size_t> index
    ) const;

    void gridSplit(
        std::vector<RadiusRegion>& regions,
        size_t k,
        size_t index = 0
    ) const;

    // ===== СТАРЫЙ API (ДЛЯ СОВМЕСТИМОСТИ) =====

    // search.cpp и codec.cpp ожидают ЭТИ имена:

    const std::vector<Interval>& GetIntervals() const {
        return getIntervals();
    }

    bool IsNarrowEnough(BaseType lower_bound) const {
        return isNarrowEnough(lower_bound);
    }

    bool IsTooWide(BaseType upper_bound) const {
        return isTooWide(upper_bound);
    }

    void Split(
        std::vector<RadiusRegion>& regions,
        size_t k,
        std::optional<size_t> index
    ) const {
        split(regions, k, index);
    }

    void GridSplit(
        std::vector<RadiusRegion>& regions,
        size_t k,
        size_t index = 0
    ) const {
        gridSplit(regions, k, index);
    }
};



struct RadiiCompare {

    bool operator()(const std::vector<Interval>& a,
                    const std::vector<Interval>& b) const;
};


/// Spiral operator lookup table

using SpiralOpRef = std::reference_wrapper<SpiralOp>;

class SpiralOpCache {

private:

    std::vector<SpiralOp> values;
    std::vector<bool> presence;

    SpiralOp identity;

    inline size_t getIndex(size_t i, size_t j, size_t k);

public:

    std::vector<Interval> radii;

    SpiralOpRef operator()(size_t base_type,
                           size_t prev_type,
                           size_t next_type);

    SpiralOpRef operator()();

    SpiralOpCache(const std::vector<Interval>& radii_);
};



struct ClockwiseDiskCompare {

    Interval center_x;
    Interval center_y;

    ClockwiseDiskCompare(const Disk& base);

    bool operator()(const DiskPtr a, const DiskPtr b) const;
};


} // namespace diskpack