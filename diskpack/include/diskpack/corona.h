#pragma once

#include <diskpack/geometry.h>

#include <queue>
#include <stack>
#include <set>
#include <list>
#include <optional>

namespace diskpack {

  class Corona;
  class CoronaSignature;
  class ConnectivityGraph;

  using CoronaSignaturePtr = std::shared_ptr<CoronaSignature>;
  using CoronaSignatureList = std::list<CoronaSignaturePtr>;

  bool CoronaSignaturePtrCompare(const CoronaSignaturePtr &a, 
                                 const CoronaSignaturePtr &b);

  bool operator< (const CoronaSignaturePtr &a, 
                  const CoronaSignaturePtr &b);   
             

  class ConnectivityGraph {     // Checks conditions required for compact packings. Used in search to validate regions.

  public:

    const size_t MAX_SIGNATURES = 5000;

    BaseType PRECISION_THRESHOLD;

    using DiffStack = std::stack<std::shared_ptr<CoronaSignatureList>>;

  private:
  
    bool broken_state = false;

    std::queue<std::tuple<size_t, size_t, size_t>> redundant_triangles;

    std::vector<CoronaSignatureList> signatures;

    std::vector<DiffStack> diffs;

    std::vector<std::vector<size_t>> transitions;

    std::vector<std::vector<bool>> edges;

    void Push(const CoronaSignature& signature);

    void Pop(const CoronaSignature& signature);

    void RemoveRedundantTriangles(std::shared_ptr<std::vector<std::shared_ptr<CoronaSignatureList>>> diff = std::shared_ptr<std::vector<std::shared_ptr<CoronaSignatureList>>>());

    void UpdateEdges();

    size_t GetTransitionConst(size_t base, 
                               size_t i, 
                               size_t j) const;

    size_t& GetTransition(size_t base, 
                           size_t i, 
                           size_t j);
                           
    void FillCorona(
      Corona& corona,
      std::list<DiskPtr> &packing, 
      size_t start_index, 
      std::set<CoronaSignaturePtr, decltype(&CoronaSignaturePtrCompare)> &unique_signatures
    );

  public:
    size_t Size() const;
    void DisplaySignatures() const;
    ConnectivityGraph(SpiralOpCache &lookup_table);
    bool HasOverflow() const;

    void Refine(SpiralOpCache &lookup_table);
    bool Restore();

    bool IsViable() const;
    bool HasTriangle(size_t i, 
                     size_t j, 
                     size_t k) const;
  };

  
  class Corona {
    const size_t DEFAULT_OPERATOR_CAPACITY = 12;

    friend class CoronaSignature;
    friend class ConnectivityGraph;

    const Disk &base;     // Central disk
    std::deque<DiskPtr> corona;     // Surrounding disks; consecutive ones are tangent


    std::vector<SpiralOpRef> operators_front;  // Auxiliary objects for computing positions of new disks
    std::vector<SpiralOpRef> operators_back;   
    IntervalPair leaf_front;              
    IntervalPair leaf_back;               
    SpiralOpCache &lookup_table;    
    std::stack<bool> push_history;       

    SpiralOp
    computeOperatorsProduct(
        const size_t &begin, 
        const size_t &end,
        const std::vector<SpiralOpRef> &operators
    ) const;

    void buildSortedCorona(     // Builds a counterclockwise ordered corona from scratch
      const std::list<DiskPtr> &packing
    );

  public:

    Corona(
      const Disk &b, 
      const std::list<DiskPtr> &packing,
          SpiralOpCache &lookup_table);

    bool isCompleted();     // Checks if corona is complete

    bool isContinuous() const;     // Checks continuity

    bool peekNewDisk(     // Computes new disk without changing state
      Disk &new_disk, 
      size_t index, 
      const std::optional<ConnectivityGraph> &graph = std::nullopt
      );       

    void Push(     // Adds disk to corona
      const DiskPtr &disk, 
      size_t index
    );

    void Pop();     // Removes last disk

    const Disk &getBase();

    void DisplaySignature();
  };


  class CoronaSignature {     // Represents adjacency between disk types in a corona. Equal signatures mean equivalent coronas.
    const size_t radii_count;
    const size_t base;
    std::vector<size_t> transitions;
    size_t& GetTransition(size_t i, 
                          size_t j);                     
    size_t GetTransitionConst(size_t i, 
                              size_t j) const;

    friend class ConnectivityGraph;
    friend Corona;

  public:
    std::vector<bool> disk_presence;
    std::vector<size_t> specimen_indices;
    CoronaSignature(Corona &specimen);
    size_t getBase() const;
    bool TestRadii(SpiralOpCache& lookup_table) const;
    bool operator<(const CoronaSignature &other) const;
    bool operator==(const CoronaSignature &other) const;
  };

}