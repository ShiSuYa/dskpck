#include <diskpack/corona.h>

#include <functional>
#include <iostream>
#include <numeric>

namespace diskpack {


SpiralOp
Corona::computeOperatorsProduct(
  const size_t &begin, 
  const size_t &end,
  const std::vector<SpiralOpRef> &operators
) const {

  switch (end - begin) {

  case 0:
    return SpiralOp();

  case 1:
    return operators[begin];

  case 2:
    return operators[begin].get() * 
           operators[begin + 1].get();

  case 3:
    return operators[begin].get() *
           (operators[begin + 1].get() * 
            operators[begin + 2].get());

  case 4:
    return (operators[begin].get() * 
            operators[begin + 1].get()) *
           (operators[begin + 2].get() * 
            operators[begin + 3].get());
  }

  auto mid = (begin + end) / 2;

  return computeOperatorsProduct(begin, mid, operators) *
         computeOperatorsProduct(mid, end, operators);
}

bool Corona::isCompleted() {

  assert(!corona.empty());

  auto cross_product =
      (corona.back()->getCenterX() - base.getCenterX()) *
      (corona.front()->getCenterY() - base.getCenterY()) -
      (corona.front()->getCenterX() - base.getCenterX()) *
      (corona.back()->getCenterY() - base.getCenterY());

  if (!(corona.size() > 2) || 
      !corona.back()->tangent(*corona.front()) ||
      !cergt(cross_product, 0.0L)) {
      return false;
  }

  bool use_front = 
      corona.front()->precision() < 
      corona.back()->precision();

  Disk new_disk;

  const Disk &old_disk = 
     *(!use_front ? corona.front() : corona.back());

  peekNewDisk(new_disk, old_disk.getType());

  bool check_intersect =
      !empty(intersect(new_disk.getCenterX(), old_disk.getCenterX())) &&
      !empty(intersect(new_disk.getCenterY(), old_disk.getCenterY()));
  return check_intersect;
}

void Corona::buildSortedCorona(const std::list<DiskPtr> &packing) {

  assert(corona.empty());

  for (auto &disk : packing) {
    if (base.tangent(*disk)) {
      corona.push_back(disk);
    }
  }

  std::sort(
    corona.begin(), 
    corona.end(), 
    ClockwiseDiskCompare{base}
  );

  for (size_t i = 0; i + 1 < corona.size(); ++i) {
    if (!corona[i]->tangent(*corona[i + 1])) {
      std::rotate(
        corona.begin(), 
        corona.begin() + i + 1, 
        corona.end()
      );
      return;
    }
  }
}

bool Corona::isContinuous() const {
  for (size_t i = 0; i + 1 < corona.size(); ++i) {
    if (!corona[i]->tangent(*corona[i + 1])) {
      return false;
    }
  }
  return true;
}

Corona::Corona(const Disk &b, const std::list<DiskPtr> &packing,
               SpiralOpCache &lookup_table_)
    : base(b), lookup_table(lookup_table_) {

  operators_back.reserve(DEFAULT_OPERATOR_CAPACITY);
  operators_front.reserve(DEFAULT_OPERATOR_CAPACITY);

  buildSortedCorona(packing);

  assert(!corona.empty());

  leaf_front =
      IntervalPair{
        corona.front()->getCenterX() - base.getCenterX(),
        corona.front()->getCenterY() - base.getCenterY()};

  leaf_back = 
      IntervalPair{
        corona.back()->getCenterX() - base.getCenterX(),
        corona.back()->getCenterY() - base.getCenterY()};
}

bool Corona::peekNewDisk(
  Disk &new_disk, 
  size_t index,
  const std::optional<ConnectivityGraph> &graph
) {

  bool use_front = 
    corona.front()->precision() < 
    corona.back()->precision();

  if (graph.has_value()) {
    if (!graph->HasTriangle(
            base.getType(),
            use_front ? corona.front()->getType() : corona.back()->getType(),
            index)) {
        return false;
    }
  }

  auto &operators = (use_front) ? operators_front : operators_back;

  auto &start_leaf = (use_front) ? leaf_front : leaf_back;

  operators.push_back(
    lookup_table(
      base.getType(),
      use_front ? corona.front()->getType() : corona.back()->getType(),
      index
    )
  );

  auto op = computeOperatorsProduct(0, operators.size(), operators);

  op.y *= (use_front ? -1 : 1);

  auto new_center = op * start_leaf;

  new_disk = Disk(
    new_center.first + base.getCenterX(),
    new_center.second + base.getCenterY(),
    lookup_table.radii[index], 
    index
  );

  operators.pop_back();

  return true;

}

void Corona::Push(const DiskPtr &disk, size_t index) {

  bool use_front = 
      corona.front()->precision() < 
      corona.back()->precision();

  push_history.push(use_front);

  auto &operators = (use_front) ? operators_front : operators_back;

  operators.push_back(
    lookup_table(
      base.getType(),
      use_front ? corona.front()->getType() : corona.back()->getType(),
      index));
  (use_front) ? corona.push_front(disk) : corona.push_back(disk);
}

void Corona::Pop() {

  auto use_front = push_history.top();
  push_history.pop();

  if (use_front) {
    operators_front.pop_back();
    corona.pop_front();
  } else {
    operators_back.pop_back();
    corona.pop_back();
  }
}

size_t ConnectivityGraph::Size() const {
  size_t result = 0;
  for (size_t base = 0; base < edges.size(); ++base) {
    result += signatures[base].size();
  }
  return result;
}

const Disk &Corona::getBase() { 
  return base; 
}

void Corona::DisplaySignature() {

  std::cerr << "signature: " << "\n";

  CoronaSignature signature(*this);

  std::cerr << "base: " << signature.base << "\n";

  for (size_t i = 0; i < lookup_table.radii.size(); ++i) {
    for (size_t j = 0; j < lookup_table.radii.size(); ++j) {
      std::cerr << signature.GetTransition(i, j) << " ";
    }
    std::cerr << "\n";
  }

  for (auto &i : signature.specimen_indices) {
    std::cerr << i << " ";
  }

  std::cerr << "\n";
}


size_t &CoronaSignature::GetTransition(size_t i, size_t j) {
  return transitions[i > j ? (i * i + i) / 2 + j : (j * j + j) / 2 + i];
}

size_t CoronaSignature::GetTransitionConst(size_t i, size_t j) const {
  return transitions[i > j ? (i * i + i) / 2 + j : (j * j + j) / 2 + i];
}

bool CoronaSignature::operator<(const CoronaSignature &other) const {
  for (size_t i = 0; i < transitions.size(); ++i) {
    if (transitions[i] != other.transitions[i]) {
      return transitions[i] < other.transitions[i];
    }
  }
  return false;
}

bool CoronaSignature::operator==(const CoronaSignature &other) const {
  for (size_t i = 0; i < transitions.size(); ++i) {
    if (transitions[i] != other.transitions[i]) {
      return false;
    }
  }
  return true;
}

CoronaSignature::CoronaSignature(Corona &specimen)
    : base(specimen.base.getType()),
      radii_count(specimen.lookup_table.radii.size()),
      specimen_indices(specimen.corona.size()),
      transitions((radii_count * radii_count + radii_count) / 2),
      disk_presence(radii_count, false) {

  if (!specimen.isCompleted()) {
    throw std::runtime_error("CoronaSignature requires completed corona");
  }

  for (size_t i = 0; i < specimen_indices.size(); ++i) {
    specimen_indices[i] = specimen.corona[i]->getType();
    disk_presence[specimen_indices[i]] = true;
  }

  for (size_t i = 0; i + 1 < specimen_indices.size(); ++i) {
    ++GetTransition(specimen_indices[i], 
                    specimen_indices[i + 1]);
  }
  ++GetTransition(specimen_indices[specimen_indices.size() - 1],
                   specimen_indices[0]);
}

bool CoronaSignature::TestRadii(SpiralOpCache &lookup_table) const {

  auto &new_radii = lookup_table.radii;

  Disk b(0, 0, new_radii[base], base);

  std::list<DiskPtr> packing;

  packing.push_back(
    std::make_shared<Disk>(
      std::move(Disk(new_radii[base] + new_radii[specimen_indices[0]], 0,
                     new_radii[specimen_indices[0]], specimen_indices[0]))));

  Corona test(b, packing, lookup_table);

  size_t cur_front = 0, cur_back = 0;

  for (size_t i = 1; i < specimen_indices.size(); ++i) {
    auto use_front =
        test.corona.front()->precision() < test.corona.back()->precision();
    if (use_front) {
      cur_front =
          (cur_front == 0 ? specimen_indices.size() - 1 : cur_front - 1);
    } else {
      cur_back = (cur_back == specimen_indices.size() - 1 ? 0 : cur_back + 1);
    }
    size_t next = specimen_indices[(use_front ? cur_front : cur_back)];
    Disk new_disk;
    test.peekNewDisk(new_disk, next);
    if (std::any_of(packing.begin(), packing.end(),
                    [&new_disk](const DiskPtr &disk) {
                      return new_disk.intersects(*disk);
                    })) {
      return false;
    }
    packing.push_back(std::make_shared<Disk>(std::move(new_disk)));
    test.Push(packing.back(), next);
  }
  return test.isCompleted();
}

size_t CoronaSignature::getBase() const { 
  return base; 
}

bool CoronaSignaturePtrCompare(const CoronaSignaturePtr &a,
                               const CoronaSignaturePtr &b) {
  return *a < *b;
}


void ConnectivityGraph::FillCorona(
    Corona &corona, 
    std::list<DiskPtr> &packing, 
    size_t start_index,
    std::set<CoronaSignaturePtr, decltype(&CoronaSignaturePtrCompare)> &unique_signatures) {
  auto &radii = corona.lookup_table.radii;
  if (corona.isCompleted()) {
    CoronaSignaturePtr signature =
        std::make_shared<CoronaSignature>(corona);
    if (unique_signatures.insert(signature).second) {
      signatures[corona.base.getType()].push_back(signature);
      Push(*signature);
    }
    return;
  }

  Disk new_disk;

  for (size_t i = start_index; i < radii.size(); ++i) {
    corona.peekNewDisk(new_disk, i);

    if (std::any_of(packing.begin(), packing.end(),
                    [&new_disk](const DiskPtr &disk) {
                      return new_disk.intersects(*disk);
                    })) {
      continue;
    }
    if (new_disk.precision() > PRECISION_THRESHOLD) {
      broken_state = true;
      return;
    }

    packing.push_back(std::make_shared<Disk>(std::move(new_disk)));
    corona.Push(packing.back(), i);

    FillCorona(corona, packing, start_index, unique_signatures);

    corona.Pop();
    packing.pop_back();
    if (HasOverflow()) {
      return;
    }
  }
}

ConnectivityGraph::ConnectivityGraph(SpiralOpCache &lookup_table)
    : diffs(lookup_table.radii.size()), 
      signatures(lookup_table.radii.size()),
      transitions(lookup_table.radii.size(),
                  std::vector<size_t>(
                      (lookup_table.radii.size() * lookup_table.radii.size() + lookup_table.radii.size()) / 2,
                      0)),
      edges(lookup_table.radii.size(),
            std::vector<bool>(lookup_table.radii.size(), false)) 
{
  std::list<DiskPtr> packing;
  auto &radii = lookup_table.radii;

  PRECISION_THRESHOLD =
      std::min_element(radii.begin(), radii.end(),
                       [](const Interval &x, const Interval &y) {
                         return x.lower() < y.lower();
                       })
          ->lower() /
      3.0L;

  for (size_t base = 0; base < radii.size(); ++base) {

    Disk base_disk(zero, zero, radii[base], base);

    std::set<CoronaSignaturePtr, decltype(&CoronaSignaturePtrCompare)>
        unique_signatures(CoronaSignaturePtrCompare);

    for (size_t start_index = 0; start_index < radii.size();
         ++start_index) {
      packing.push_back(
          std::make_shared<Disk>(
            radii[base] + radii[start_index], zero,
            radii[start_index], start_index));

      Corona corona(base_disk, packing, lookup_table);

      FillCorona(corona, packing, start_index, unique_signatures);

      packing.clear();
    }
  }

  for (size_t base = 0; base < edges.size(); ++base) {
    for (size_t i = 0; i < edges.size(); ++i) {
      for (size_t j = i; j < edges.size(); ++j) {
        if (GetTransitionConst(base, i, j) == 0) {
          redundant_triangles.push(std::make_tuple(base, i, j));
        }
      }
    }
  }

  if (HasOverflow()) {
    return;
  }

  RemoveRedundantTriangles();
  UpdateEdges();
}

void ConnectivityGraph::DisplaySignatures() const {
  for (size_t base = 0; base < edges.size(); ++base) {
    for (auto signature : signatures[base]) {
      std::cerr << base << " : ";
      for (auto index : signature->specimen_indices) {
        std::cerr << index << " ";
      }
      std::cerr << "\n";
    }
  }
}

void ConnectivityGraph::Refine(SpiralOpCache &lookup_table) {
  auto &radii = lookup_table.radii;

  std::shared_ptr<std::vector<std::shared_ptr<CoronaSignatureList>>>
      invalid_signatures =
          std::make_shared<std::vector<std::shared_ptr<CoronaSignatureList>>>(
              diffs.size());
  for (auto &list : *invalid_signatures) {
    list = std::make_shared<CoronaSignatureList>();
  }

  for (size_t base = 0; base < radii.size(); ++base) {
    for (auto it = signatures[base].begin(); it != signatures[base].end();) {
      auto signature = *it;
      if (signature->TestRadii(lookup_table)) {
        ++it;
        continue;
      }
      invalid_signatures->operator[](base)->push_back(signature);
      it = signatures[base].erase(it);
      Pop(*signature);
    }
  }
  RemoveRedundantTriangles(invalid_signatures);
  for (size_t base = 0; base < edges.size(); ++base) {
    diffs[base].push(invalid_signatures->operator[](base));
  }
  UpdateEdges();
}

bool ConnectivityGraph::Restore() {
  if (std::any_of(diffs.begin(), diffs.end(),
                  [](const DiffStack &s) { return s.empty(); })) {
    return false;
  }

  for (size_t base = 0; base < diffs.size(); ++base) {
    auto diff = diffs[base].top();
    diffs[base].pop();
    for (auto signature : *diff) {
      Push(*signature);
    }
    signatures[base].splice(signatures[base].end(), std::move(*diff));
  }

  UpdateEdges();
  return true;
}

bool ConnectivityGraph::HasOverflow() const {
  return broken_state || std::any_of(signatures.begin(), signatures.end(),
                                  [this](const CoronaSignatureList &x) {
                                    return x.size() > MAX_SIGNATURES;
                                  });
}

void ConnectivityGraph::Push(const CoronaSignature &signature) {
  for (size_t x = 0; x < signature.transitions.size(); ++x) {
    transitions[signature.base][x] += signature.transitions[x];
  }
}

void ConnectivityGraph::Pop(const CoronaSignature &signature) {
  for (size_t i = 0; i < edges.size(); ++i) {
    for (size_t j = i; j < edges.size(); ++j) {
      GetTransition(signature.base, i, j) -=
          signature.GetTransitionConst(i, j);
      if (signature.GetTransitionConst(i, j) > 0 &&
          GetTransitionConst(signature.base, i, j) == 0) {
        redundant_triangles.push(std::make_tuple(i, j, signature.base));
      }
    }
  }
}

void ConnectivityGraph::RemoveRedundantTriangles(
    std::shared_ptr<std::vector<std::shared_ptr<CoronaSignatureList>>> diff) 
{
  std::vector<CoronaSignatureList> removed_signatures(edges.size());

  while (!redundant_triangles.empty()) {

    auto [i, j, k] = redundant_triangles.front();

    redundant_triangles.pop();

    std::vector<size_t> ijk{i, j, k};

    for (size_t rotation_index = 0; rotation_index < ijk.size();
         ++rotation_index) {
      std::rotate(ijk.begin(), std::next(ijk.begin()), ijk.end());
      if (GetTransitionConst(ijk[0], ijk[1], ijk[2]) == 0) {
        continue;
      }
      for (auto it = signatures[ijk[0]].begin();
           it != signatures[ijk[0]].end();) {
        auto signature = *it;
        if (signature->GetTransitionConst(ijk[1], ijk[2]) > 0) {
          removed_signatures[ijk[0]].push_back(signature);
          it = signatures[ijk[0]].erase(it);
          Pop(*signature);
        } else {
          ++it;
        }
      }
    }
  }

  if (diff) {
    for (size_t base = 0; base < edges.size(); ++base) {
      diff->operator[](base)->splice(diff->operator[](base)->end(),
                                     std::move(removed_signatures[base]));
    }
  }
}

size_t &ConnectivityGraph::GetTransition(size_t base, size_t i, size_t j) {
  return transitions[base][i > j ? (i * i + i) / 2 + j : (j * j + j) / 2 + i];
}

size_t ConnectivityGraph::GetTransitionConst(size_t base, size_t i,
                                              size_t j) const {
  return transitions[base][i > j ? (i * i + i) / 2 + j : (j * j + j) / 2 + i];
}

void ConnectivityGraph::UpdateEdges() {

  for (size_t i = 0; i < edges.size(); ++i) {
    for (size_t j = i; j < edges.size(); ++j) {

      edges[i][j] = false;

      for (size_t x = 0; x < edges.size(); ++x) {
        
        edges[i][j] = edges[i][j] || (GetTransitionConst(i, j, x) > 0 &&
                                      GetTransitionConst(j, i, x) > 0);
        edges[j][i] = edges[i][j];
      }
    }
  }
}

bool ConnectivityGraph::IsViable() const {
  // Build a graph of disk types.
  // An edge means that the corresponding disk types can appear together in at least one locally valid corona configuration.

  // Condition 1: the graph of disk types must be connected
  
  {
    std::vector<bool> unvisited(edges.size(), true);
    std::queue<size_t> q;

    q.push(0);
    unvisited[0] = false;

    while (!q.empty()) {
      auto vertice = q.front();
      q.pop();

      for (size_t neighbor = 0; neighbor < edges.size(); ++neighbor) {
        if (unvisited[neighbor] && edges[vertice][neighbor]) {
          unvisited[neighbor] = false;
          q.push(neighbor);
        }
      }
    }

    if (std::any_of(unvisited.begin(), unvisited.end(),
                    [](bool b) { return b; })) {
      return false;
    }
  }

  /// Condition 2: if one disk type is removed, its coronas must connect all remaining components

  {
    for (size_t removed = 0; removed < edges.size();
         ++removed) {
      std::vector<size_t> components(edges.size(), 0);
      components[removed] = edges.size();

      size_t current_component = 0;

      for (size_t vertice = 0; vertice < edges.size(); ++vertice) {
        if (components[vertice]) {
          continue;
        }

        ++current_component;

        std::queue<size_t> q;
        q.push(vertice);
        components[vertice] = current_component;

        while (!q.empty()) {
          size_t x = q.front();
          q.pop();

          for (size_t neighbor = 0; neighbor < edges.size(); ++neighbor) {
            if (components[neighbor] || !edges[x][neighbor]) {
              continue;
            }
            components[neighbor] = current_component;
            q.push(neighbor);
          }
        }
      }

      if (!std::any_of(signatures[removed].begin(),
                       signatures[removed].end(),
                       [current_component, this,
                        &components](const CoronaSignaturePtr &signature) {
                         for (size_t component = 1; component <= current_component;
                              ++component) {
                           bool connects = false;
                           for (size_t vertice = 0; vertice < edges.size();
                                ++vertice) {
                             connects =
                                 connects ||
                                 (signature->disk_presence[vertice] &&
                                  components[vertice] == component);
                           }
                           if (!connects) {
                             return false;
                           }
                         }
                         return true;
                       })) {
        return false;
      }
    }
  }

  // Condition 3: compatible coronas must connect all disk types

  {
    std::vector<size_t> components(edges.size() * edges.size() * edges.size());
    std::vector<size_t> component_size(components.size(), 1);
    std::iota(components.begin(), components.end(), 0);

    auto getIndex = [&](size_t i, size_t j, size_t k) {
      std::vector<size_t> ijk{i, j, k};
      std::sort(ijk.begin(), ijk.end());
      return ijk[0] + edges.size() * ijk[1] +
             (edges.size() * edges.size()) * ijk[2];
    };

    std::function<size_t(size_t)> find = [&](size_t x) {
      components[x] = (components[x] == x ? x : find(components[x]));
      return components[x];
    };

    std::function<void(size_t, size_t size_t)> unite = [&](size_t x, size_t y) {
      x = find(x);
      y = find(y);
      if (x != y) {
        if (component_size[x] > component_size[y]) {
          std::swap(x, y);
        }
        components[x] = y;
        component_size[y] += component_size[x];
      }
    };

    for (size_t base = 0; base < edges.size(); ++base) {
      for (auto signature : signatures[base]) {
        size_t first_i = signature->specimen_indices[0],
               first_j = signature->specimen_indices[1];
        size_t first_index = getIndex(base, first_i, first_j);

        for (size_t i = 0; i < edges.size(); ++i) {
          for (size_t j = i; j < edges.size(); ++j) {
            if (signature->GetTransitionConst(i, j)) {
              unite(first_index, getIndex(base, i, j));
            }
          }
        }
      }
    }

    std::vector<std::vector<bool>> cdt(
        edges.size() * edges.size() * edges.size(),
        std::vector<bool>(edges.size(), false));
    for (size_t i = 0; i < edges.size(); ++i) {
      for (size_t j = i; j < edges.size(); ++j) {
        for (size_t k = j; k < edges.size(); ++k) {
          size_t component = find(getIndex(i, j, k));

          cdt[component][i] = true;
          cdt[component][j] = true;
          cdt[component][k] = true;
        }
      }
    }

    if (!std::any_of(cdt.begin(), cdt.end(),
                     [](const std::vector<bool> &disks) {
                       return !std::any_of(disks.begin(), disks.end(),
                                           [](bool u) { return !u; });
                     })) {
      return false;
    }
  }

  return true;
}

bool ConnectivityGraph::HasTriangle(size_t i, size_t j, size_t k) const {
  return GetTransitionConst(i, j, k) && 
         GetTransitionConst(k, i, j) &&
         GetTransitionConst(j, k, i);
}

}