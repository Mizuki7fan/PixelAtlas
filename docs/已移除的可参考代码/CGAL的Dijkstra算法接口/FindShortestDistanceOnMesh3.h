#pragma once
#include <CGAL-Interface/CGAL-Interface.h>
#include <boost/graph/adjacency_list.hpp>
#include <boost/graph/dijkstra_shortest_paths.hpp>

#include <boost/range/any_range.hpp>
template <typename Type, typename Traversal>
struct AnyRange : boost::any_range<Type, Traversal> {
  typedef boost::any_range<Type, Traversal> BaseType;

  template <typename IteratorType> struct IteratorTypeChecker {
    typedef
        typename std::iterator_traits<IteratorType>::reference ReferenceType;
    typedef typename std::iterator_traits<IteratorType>::value_type ValueType;

    constexpr static bool value = !std::is_same_v<ReferenceType, ValueType>;
  };

  template <class WrappedRange>
  AnyRange(WrappedRange &wrapped_range) : BaseType(wrapped_range) { // NOLINT
    static_assert(
        IteratorTypeChecker<decltype(std::ranges::begin(wrapped_range))>::value,
        "This range is not allowed to be converted into AnyRange.");
  }

  template <class WrappedRange>
  AnyRange(const WrappedRange &wrapped_range)
      : BaseType(wrapped_range) { // NOLINT
    static_assert(
        IteratorTypeChecker<decltype(std::ranges::begin(wrapped_range))>::value,
        "This range is not allowed to be converted into AnyRange.");
  }

  template <class Iterator>
  AnyRange(Iterator first, Iterator last) : BaseType(first, last) {
    static_assert(IteratorTypeChecker<Iterator>::value,
                  "This range is not allowed to be converted into AnyRange.");
  }

  AnyRange() : BaseType() {}
};
template <typename Type>
using AnySinglePassRange = AnyRange<Type, boost::single_pass_traversal_tag>;

class DijkstraAllTargetReachGuardingVisitor : public boost::dijkstra_visitor<> {
public:
  explicit DijkstraAllTargetReachGuardingVisitor(
      const AnySinglePassRange<CGAL::SM_Vertex_index> &target_vtxs);

  void examine_vertex(CGAL::SM_Vertex_index vertex, const cgl::SurfaceMesh3 &);

  // This exception is used to stop the Dijkstra algorithm when all targets are
  // reached.
  struct AllTargetsReachedException {};

protected:
  std::set<CGAL::SM_Vertex_index> unreached_target_vtxs_;
};

using PredecessorMap = boost::iterator_property_map<
    std::vector<CGAL::SM_Vertex_index>::iterator,
    boost::typed_identity_property_map<CGAL::SM_Vertex_index>>;

using DistanceMap = boost::iterator_property_map<
    std::vector<double>::iterator,
    boost::typed_identity_property_map<CGAL::SM_Vertex_index>>;

void FindShortestPathsOnMesh3(
    const cgl::SurfaceMesh3 &mesh, CGAL::SM_Vertex_index source,
    const AnySinglePassRange<CGAL::SM_Vertex_index> &targets,
    std::vector<CGAL::SM_Vertex_index> &map_vtx_to_prev_vtx,
    std::vector<double> &map_vtx_to_distance);

std::pair<double, std::vector<CGAL::SM_Halfedge_index>>
FindShortestHalfedgePathOnMesh3(const cgl::SurfaceMesh3 &mesh,
                                cgl::SurfaceMesh3::vertex_index source,
                                cgl::SurfaceMesh3::vertex_index target);