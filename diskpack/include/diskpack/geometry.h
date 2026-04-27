#pragma once

#include <boost/numeric/interval.hpp>
#include <memory>
#include <optional>
#include <vector>

namespace diskpack {

using namespace boost::numeric;

using BaseType = long double;     // Numeric type alias, can be replaced with higher precision type

using Interval = interval<     // Interval arithmetic   
    BaseType,
    interval_lib::policies<
        interval_lib::save_state<interval_lib::rounded_transc_exact<BaseType>>,
        interval_lib::checking_base<BaseType>>>;

using IntervalPair = std::pair<Interval, Interval>;

const BaseType epsilon = .00000000001;

const Interval one{1 - epsilon, 1 + epsilon}, zero{-epsilon, epsilon};

struct SpiralOp {     // Computes the position of a disk tangent to two given disks

  Interval x;
  Interval y;

  SpiralOp(const Interval &x_, const Interval &y_);
  SpiralOp(const SpiralOp &other);
  SpiralOp(SpiralOp &&other);

  SpiralOp &operator=(const SpiralOp &other);
  SpiralOp &operator=(SpiralOp &&other);

  SpiralOp(const Interval &base_r, const Interval &prev_r,
           const Interval &next_r, const size_t &base_t,
           const size_t &prev_t, const size_t &next_t);

  IntervalPair operator*(const IntervalPair &vec) const noexcept;
  SpiralOp
  operator*(const SpiralOp &other) const noexcept;
  SpiralOp();
};

class Disk {
  Interval center_x, center_y, radius;
  size_t disk_type;

  public:
  Disk(Interval center_x_, Interval center_y_, Interval radius_,
       size_t disk_type_);
  Disk();
  Disk(const Disk &other) noexcept;
  Disk(Disk &&other) noexcept;
  Disk &operator=(const Disk &other);
  Disk &operator=(Disk &&other);

  Interval getNorm() const;

  const Interval &getRadius() const;
  const Interval &getCenterX() const;
  const Interval &getCenterY() const;

  size_t getType() const;

  bool intersects(const Disk &other) const;
  bool tangent(const Disk &other) const;
  bool disjoint(const Disk &other) const;

  BaseType precision() const;
};

using DiskPtr = std::shared_ptr<Disk>;

bool compareDiskNorm(const DiskPtr a, const DiskPtr b);


class RadiusRegion {      // Represents a rectangular region of radii with basic operations

  std::vector<Interval> intervals;

public:
    
    RadiusRegion(const std::vector<Interval> &intervals_);
    RadiusRegion(std::vector<Interval> &&intervals_);
    
    const std::vector<Interval>& getIntervals() const;
    
    bool isNarrowEnough(BaseType lower_bound) const;
    bool isTooWide(BaseType upper_bound) const;

    Interval getMinInterval() const;
    Interval getMaxInterval() const;

    void split(     // Divides the region into k equal parts along one coordinate
        std::vector<RadiusRegion> &regions, 
        size_t k, 
        std::optional<size_t> index
    ) const;   
                                                                                                  
    void gridSplit(     // Divides the region into k^n parts across all coordinates
        std::vector<RadiusRegion> &regions, 
        size_t k, 
        size_t index = 0
    ) const;                                                                                                         
};

struct RadiiCompare {
  bool operator()(const std::vector<Interval> &a, 
                  const std::vector<Interval> &b) const;
};


// Cache of spiral operators to avoid repeated computations

using SpiralOpRef = std::reference_wrapper<SpiralOp>;

class SpiralOpCache {

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

  SpiralOpCache(const std::vector<Interval> &radii_);
};

struct ClockwiseDiskCompare {

  Interval center_x;
  Interval center_y;

  ClockwiseDiskCompare(const Disk &base);

  bool operator()(const DiskPtr a, const DiskPtr b) const;
};

}