#include <diskpack/generator.h>

#include <iostream>
#include <algorithm>
#include <numeric>

std::random_device rd;

namespace diskpack {

const Interval zero = Interval(0.0);

BasicGenerator::BasicGenerator(
        const std::vector<Interval>& radii_,
        const RadiusType& packing_radius_,
        const RadiusType& precision_upper_bound_,
        const size_t& size_upper_bound_,
        const size_t& max_ignored_radii_)
    :
      radii(radii_),
      precision_upper_bound(precision_upper_bound_),
      packing_radius(packing_radius_),
      size_upper_bound(size_upper_bound_),
      max_ignored_radii(max_ignored_radii_),
      disk_queue(LessNormCompare()),
      op_cache(radii_),
      frequency_table(radii_.size(), 0),
      generated_radius(0),
      rng(rd()),
      graph(op_cache)
{
    if (graph->HasOverflow()) {
        graph.reset();
    }
}






PackingStatus BasicGenerator::Generate(const size_t& central_disk_type)
{
    if (radii.empty()) {
        throw std::runtime_error(
            "BasicGenerator::Generate() called with no radii");
    }

    if (graph.has_value()) {
        if (!graph->IsViable()) {
            std::cerr << "Not viable!\n";
            return PackingStatus::invalid;
        }
    }

    Reset();

    Push(
        Disk(zero, zero, radii[central_disk_type], central_disk_type),
        central_disk_type
    );

    for (size_t i = 0; i < radii.size(); ++i) {

        if (i == central_disk_type && radii.size() > 1) {
            continue;
        }

        Push(
            Disk(radii[central_disk_type] + radii[i], zero, radii[i], i),
            i
        );

        auto status = AdvancePacking();

        if (status != PackingStatus::invalid) {
            return status;
        }

        Pop(i);
    }

    Pop(central_disk_type);

    return PackingStatus::invalid;
}






PackingStatus BasicGenerator::Resume()
{
    if (packing.size() < 2) {
        return PackingStatus::invalid;
    }

    return AdvancePacking();
}






void BasicGenerator::Reset()
{
    disk_queue.clear();
    packing.clear();
    generated_radius = 0;

    std::fill(
        frequency_table.begin(),
        frequency_table.end(),
        0
    );
}






void BasicGenerator::ShuffleIndexes(std::vector<size_t>& shuffle)
{
    std::iota(shuffle.begin(), shuffle.end(), 0);

    std::shuffle(
        shuffle.begin(),
        shuffle.end(),
        rng
    );
}






/// ------------------------------------------------------------
/// Recursive search
/// ------------------------------------------------------------

PackingStatus BasicGenerator::AdvancePacking()
{
    if (disk_queue.empty()) {
        return PackingStatus::invalid;
    }

    auto base = disk_queue.extract(disk_queue.begin()).value();

    if (base->precision() > precision_upper_bound) {

        SetGeneratedRadius(*base);

        return PackingStatus::precision_error;
    }

    if (!IsInBounds(*base) || PackingIsLargeEnough()) {

        disk_queue.insert(base);

        if (!PackingSatisfiesConstraints()) {
            return PackingStatus::invalid;
        }

        SetGeneratedRadius(*base);

        return PackingStatus::complete;
    }

    Corona corona(*base, packing, op_cache);

    if (!corona.isContinuous()) {

        SetGeneratedRadius(*base);

        return PackingStatus::corona_error;
    }

    auto status = GapFill(corona);

    if (status == PackingStatus::invalid) {
        disk_queue.insert(base);
    }

    return status;
}






PackingStatus BasicGenerator::GapFill(Corona& corona)
{
    if (corona.isCompleted()) {
        return AdvancePacking();
    }

    std::vector<size_t> shuffle(radii.size());

    ShuffleIndexes(shuffle);

    Disk new_disk;

    for (size_t i = 0; i < radii.size(); ++i) {

        if (!corona.peekNewDisk(new_disk, shuffle[i], graph)) {
            continue;
        }

        if (HasIntersection(new_disk)) {
            continue;
        }

        Push(std::move(new_disk), shuffle[i]);

        corona.PushDisk(packing.back(), shuffle[i]);

        auto status = GapFill(corona);

        if (status != PackingStatus::invalid) {
            return status;
        }

        corona.PopDisk();

        Pop(shuffle[i]);
    }

    return PackingStatus::invalid;
}






std::ostream& operator<<(std::ostream& out, PackingStatus status)
{
    switch (status) {

        case PackingStatus::complete:
            return out << "complete";

        case PackingStatus::invalid:
            return out << "invalid";

        case PackingStatus::corona_error:
            return out << "corona_error";

        case PackingStatus::precision_error:
            return out << "precision_error";
    }

    return out;
}






/// ------------------------------------------------------------
/// Packing stack operations
/// ------------------------------------------------------------

void BasicGenerator::Push(Disk&& new_disk, size_t index)
{
    packing.push_back(
        std::make_shared<Disk>(std::move(new_disk))
    );

    disk_queue.insert(packing.back());

    ++frequency_table[index];
}






void BasicGenerator::Pop(size_t index)
{
    --frequency_table[index];

    disk_queue.erase(packing.back());

    packing.pop_back();
}






/// ------------------------------------------------------------
/// Geometry checks
/// ------------------------------------------------------------

bool BasicGenerator::HasIntersection(const Disk& new_disk) const
{
    return std::any_of(
        packing.begin(),
        packing.end(),
        [&new_disk](const DiskPtr& disk)
        {
            return new_disk.intersects(*disk);
        }
    );
}






bool BasicGenerator::IsInBounds(const Disk& disk) const
{
    return cerle(
        disk.getNorm(),
        (packing_radius - disk.getRadius()) *
        (packing_radius - disk.getRadius())
    );
}






bool BasicGenerator::PackingSatisfiesConstraints() const
{
    return std::count(
        frequency_table.begin(),
        frequency_table.end(),
        0
    ) <= max_ignored_radii;
}






bool BasicGenerator::PackingIsLargeEnough() const
{
    return packing.size() >= size_upper_bound;
}






/// ------------------------------------------------------------
/// Parameter setters
/// ------------------------------------------------------------

void BasicGenerator::SetGeneratedRadius(const Disk& furthest_disk)
{
    generated_radius =
        median(
            sqrt(furthest_disk.getNorm()) +
            furthest_disk.getRadius()
        );
}






void BasicGenerator::SetRadii(const std::vector<Interval>& radii_)
{
    Reset();

    radii = radii_;

    op_cache = std::move(SpiralOpCache(radii_));

    frequency_table.resize(radii.size(), 0);
}






void BasicGenerator::SetPackingRadius(const RadiusType& new_packing_radius)
{
    packing_radius = new_packing_radius;
}






void BasicGenerator::SetSizeUpperBound(const size_t& new_size)
{
    size_upper_bound = new_size;
}






/// ------------------------------------------------------------
/// Getters
/// ------------------------------------------------------------

const std::list<DiskPtr>& BasicGenerator::GetPacking()
{
    return packing;
}






const RadiusType& BasicGenerator::GetRadius()
{
    return packing_radius;
}






const RadiusType& BasicGenerator::GetGeneratedRadius()
{
    return generated_radius;
}

} // namespace diskpack