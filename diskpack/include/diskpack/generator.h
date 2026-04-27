#pragma once

#include <diskpack/corona.h>
#include <set>
#include <random>

namespace diskpack {

enum PackingStatus {
  complete,     // Packing was generated successfully
  invalid,     // No valid packing exists for these radii
  precision_error,     // Coordinate intervals became too wide
  corona_error     // Corona construction failed
};

using QueueType = std::multiset<DiskPtr, decltype(&compareDiskNorm)>;     // Ordered disk queue; returns disks closer to (0, 0) first.

std::ostream &operator<<(std::ostream &out, PackingStatus status);


class BasicGenerator {

protected:

  std::mt19937 g;

  std::vector<Interval> radii;     // Disk radii

  const BaseType precision_upper_bound;     // Maximum allowed coordinate interval width

  BaseType packing_radius;     // Target radius of the generated packing

  size_t size_upper_bound;     // Maximum number of disks

  std::list<DiskPtr> packing;     // Current packing

  QueueType disk_queue;     // Queue used during generation

  SpiralOpCache lookup_table;     // Cached geometric operators

  std::vector<size_t> frequency_table;     // Counts disks of each type

  size_t max_ignored_radii;

  BaseType generated_radius;     // Radius actually covered by the packing

  PackingStatus GapFill(Corona &corona);     // Fills gaps until the corona is complete

  PackingStatus AdvancePacking();     // Picks the next disk and fills its corona

  void ShuffleIndexes(std::vector<size_t> &shuffle);     // Randomizes disk type order

  void Push(Disk &&new_disk, size_t index);     // Adds a disk and updates generator state

  void Pop(size_t index);     // Removes the last added disk and restores state

  void SetGeneratedRadius(const Disk &furthest_disk);



  bool HasIntersection(const Disk &new_disk) const;      // Checks whether a new disk intersects the packing

  bool IsInBounds(const Disk &disk) const;     // Checks whether a disk is inside the target region

  bool PackingSatisfiesConstraints() const;     // Checks packing constraints

  bool PackingIsLargeEnough() const;      // Checks if the size limit was reached

  
public:

  std::optional<ConnectivityGraph> graph;

  BasicGenerator(  
    const std::vector<Interval> &radii,
    const BaseType &packing_radius,
    const BaseType &precision_upper_bound,
    const size_t &size_upper_bound,
    const size_t &max_ignored_radii = 0
  );

  PackingStatus Generate(const size_t &central_disk_type);     // Starts packing generation

  void Reset();     // Clears the current state

  PackingStatus Resume();     // Continues generation after parameter changes

  
  

  void SetPackingRadius(const BaseType &new_packing_radius);     // Updates target packing radius

  void SetSizeUpperBound(const size_t &new_size);     // Updates disk count limit

  void SetRadii(const std::vector<Interval> &radii_);




  const BaseType &GetGeneratedRadius();     // Returns covered packing radius

  const BaseType &GetRadius();     // Returns target packing radius

  const std::list<DiskPtr> &GetPacking();     // Returns generated packing

};

}