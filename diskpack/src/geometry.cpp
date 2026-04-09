#include <diskpack/geometry.h>
#include <algorithm>
#include <cmath>

namespace diskpack {

/// =========================
/// RadiusRegion
/// =========================

void RadiusRegion::split(std::vector<RadiusRegion>& regions,
                         size_t k,
                         std::optional<size_t> index) const {

    regions.clear();
    regions.reserve(k);

    auto intervals_copy = intervals;

    auto it = !index.has_value()
        ? std::max_element(
              intervals_copy.begin(),
              intervals_copy.end(),
              [](const Interval& a, const Interval& b) {
                  return width(a) < width(b);
              })
        : intervals_copy.begin() + index.value();

    auto initial = *it;
    auto sub_width = (initial.upper() - initial.lower()) / k;

    auto next = std::next(it);

    for (size_t i = 0; i < k; ++i) {

        auto l = initial.lower() + i * sub_width;

        auto u = (i < k - 1)
                     ? initial.lower() + (i + 1) * sub_width
                     : initial.upper();

        *it = Interval{l, u};

        if (next != intervals_copy.end()) {
            if (cergt(*it, *next))
                continue;
        }

        if (it != intervals_copy.begin()) {
            if (cergt(*std::prev(it), *it))
                continue;
        }

        regions.emplace_back(intervals_copy);
    }
}

void RadiusRegion::gridSplit(std::vector<RadiusRegion>& regions,
                             size_t k,
                             size_t index) const {

    if (index == 0) {

        regions.clear();

        size_t s = 1;

        for (size_t i = 0; i + 1 < intervals.size(); ++i)
            s *= k;

        regions.reserve(s);
    }

    if (index + 1 >= intervals.size()) {
        regions.emplace_back(intervals);
        return;
    }

    std::vector<RadiusRegion> children;

    split(children, k, {});

    for (auto& r : children)
        r.gridSplit(regions, k, index + 1);
}

const std::vector<Interval>& RadiusRegion::getIntervals() const {
    return intervals;
}

bool RadiusRegion::isNarrowEnough(BaseType lower_bound) const {

    return !std::any_of(
        intervals.begin(),
        intervals.end(),
        [&](const Interval& x) { return width(x) > lower_bound; });
}

bool RadiusRegion::isTooWide(BaseType upper_bound) const {

    return std::any_of(
        intervals.begin(),
        intervals.end(),
        [&](const Interval& x) { return width(x) > upper_bound; });
}

RadiusRegion::RadiusRegion(const std::vector<Interval>& intervals_)
    : intervals(intervals_) {}

RadiusRegion::RadiusRegion(std::vector<Interval>&& intervals_)
    : intervals(std::move(intervals_)) {}

Interval RadiusRegion::getMinInterval() const {
    return intervals.front();
}

Interval RadiusRegion::getMaxInterval() const {
    return intervals.back();
}

/// =========================
/// RadiiCompare
/// =========================

bool RadiiCompare::operator()(const std::vector<Interval>& a,
                              const std::vector<Interval>& b) const {

    for (size_t i = 0; i < a.size(); ++i) {

        if (a[i].lower() != b[i].lower())
            return a[i].lower() < b[i].lower();

        if (a[i].upper() != b[i].upper())
            return a[i].upper() < b[i].upper();
    }

    return false;
}

/// =========================
/// SpiralOp
/// =========================

SpiralOp::SpiralOp(const Interval& x_, const Interval& y_)
    : x{x_}, y{y_} {}

SpiralOp::SpiralOp(const SpiralOp& other)
    : x{other.x}, y{other.y} {}

SpiralOp::SpiralOp(SpiralOp&& other)
    : x{std::move(other.x)}, y{std::move(other.y)} {}

SpiralOp& SpiralOp::operator=(const SpiralOp& other) {

    x = other.x;
    y = other.y;

    return *this;
}

SpiralOp& SpiralOp::operator=(SpiralOp&& other) {

    x = std::move(other.x);
    y = std::move(other.y);

    return *this;
}

SpiralOp::SpiralOp(const Interval& base_r,
                   const Interval& prev_r,
                   const Interval& next_r,
                   const size_t& base_t,
                   const size_t& prev_t,
                   const size_t& next_t) {

    const auto& b_r = base_r;
    const auto& p_r = prev_r;
    const auto& n_r = next_r;

    const auto& b_t = base_t;
    const auto& p_t = prev_t;
    const auto& n_t = next_t;

    x = (p_t == b_t
             ? ONE / 2.0L
             : (b_t == n_t
                    ? ONE * 2.0L / square(ONE + p_r / b_r)
                    : (n_t == p_t
                           ? ONE - ONE * 2.0L / square(b_r / n_r + ONE)
                           : ONE / (ONE + p_r / b_r) +
                                 n_r * (ONE - p_r / b_r) /
                                     (b_r * square(ONE + p_r / b_r)))));

    auto t = n_r / (b_r + p_r);

    y = (b_t == p_t && n_t == p_t
             ? ONE * std::sqrt(3.0L) / 2.0L
             : (b_t == p_t
                    ? sqrt(t * (ONE + t))
                    : 2.0L *
                          sqrt(t * (ONE + t) /
                               (ONE * 2.0L + p_r / b_r + b_r / p_r))));
}

IntervalPair SpiralOp::operator*(const IntervalPair& vec) const noexcept {

    return {
        vec.first * x - vec.second * y,
        vec.second * x + vec.first * y};
}

SpiralOp SpiralOp::operator*(const SpiralOp& other) const noexcept {

    return SpiralOp{
        other.x * x - other.y * y,
        other.y * x + other.x * y};
}

SpiralOp::SpiralOp()
    : SpiralOp{1, 0} {}

/// =========================
/// Disk
/// =========================

Disk::Disk(Interval cx,
           Interval cy,
           Interval r,
           size_t type)
    : center_x{cx},
      center_y{cy},
      radius{r},
      disk_type{type} {}

Disk::Disk(const Disk& other) noexcept
    : center_x{other.center_x},
      center_y{other.center_y},
      radius{other.radius},
      disk_type{other.disk_type} {}

Disk::Disk(Disk&& other) noexcept
    : center_x{std::move(other.center_x)},
      center_y{std::move(other.center_y)},
      radius{std::move(other.radius)},
      disk_type{std::move(other.disk_type)} {}

Disk& Disk::operator=(const Disk& other) {

    center_x = other.center_x;
    center_y = other.center_y;
    radius = other.radius;
    disk_type = other.disk_type;

    return *this;
}

Disk& Disk::operator=(Disk&& other) {

    center_x = std::move(other.center_x);
    center_y = std::move(other.center_y);
    radius = std::move(other.radius);
    disk_type = std::move(other.disk_type);

    return *this;
}

Disk::Disk()
    : center_x(0, 0),
      center_y(0, 0),
      radius(0, 0),
      disk_type(0) {}

Interval Disk::getNorm() const {
    return square(center_x) + square(center_y);
}

const Interval& Disk::getRadius() const {
    return radius;
}

const Interval& Disk::getCenterX() const {
    return center_x;
}

const Interval& Disk::getCenterY() const {
    return center_y;
}

size_t Disk::getType() const {
    return disk_type;
}

BaseType Disk::precision() const {
    return std::min(width(center_x), width(center_y));
}

/// =========================
/// Disk geometry
/// =========================

inline Interval diskGap(const Disk& a, const Disk& b) {

    return sqrt(square(a.getCenterX() - b.getCenterX()) +
                square(a.getCenterY() - b.getCenterY())) -
           a.getRadius() - b.getRadius();
}

bool Disk::intersects(const Disk& other) const {
    return cerlt(diskGap(*this, other), 0.0L);
}

bool Disk::tangent(const Disk& other) const {
    return zero_in(diskGap(*this, other));
}

bool Disk::disjoint(const Disk& other) const {
    return cergt(diskGap(*this, other), 0.0L);
}

/// =========================
/// Comparators
/// =========================

bool compareDiskNorm(const DiskPtr a, const DiskPtr b) {

    return median(sqrt(a->getNorm()) + a->getRadius()) <
           median(sqrt(b->getNorm()) + b->getRadius());
}

ClockwiseDiskCompare::ClockwiseDiskCompare(const Disk& base) {

    center_x = base.getCenterX();
    center_y = base.getCenterY();
}

bool ClockwiseDiskCompare::operator()(const DiskPtr a,
                                      const DiskPtr b) const {

    auto ax = a->getCenterX() - center_x;
    auto ay = a->getCenterY() - center_y;

    auto bx = b->getCenterX() - center_x;
    auto by = b->getCenterY() - center_y;

    return (
        cergt(ax, 0.0L) && zero_in(ay)
            ? true
            : (cergt(bx, 0.0L) && zero_in(by)
                   ? false
                   : (cerlt(ay * by, 0.0L)
                          ? cergt(ay, 0.0L)
                          : cergt(by * ax - ay * bx, 0.0L))));
}

/// =========================
/// SpiralOpCache
/// =========================

inline size_t SpiralOpCache::getIndex(size_t i,
                                      size_t j,
                                      size_t k) {

    return i +
           j * radii.size() +
           k * radii.size() * radii.size();
}

SpiralOpRef SpiralOpCache::operator()(size_t base_t,
                                      size_t prev_t,
                                      size_t next_t) {

    auto index = getIndex(base_t, prev_t, next_t);

    if (!presence[index]) {

        values[index] = SpiralOp{
            radii[base_t],
            radii[prev_t],
            radii[next_t],
            base_t,
            prev_t,
            next_t};

        presence[index] = true;
    }

    return std::ref(values[index]);
}

SpiralOpCache::SpiralOpCache(const std::vector<Interval>& radii_)
    : radii(radii_),
      identity(),
      values(radii_.size() * radii_.size() * radii_.size()),
      presence(radii_.size() * radii_.size() * radii_.size()) {}

SpiralOpRef SpiralOpCache::operator()() {
    return std::ref(identity);
}

} 