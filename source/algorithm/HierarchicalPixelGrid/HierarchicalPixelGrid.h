#pragma once
#include <array>
#include <fstream>
#include <vector>

class GridFace;
class GridEdge;
class GridHalfEdge;
class GridVertex;

class GridFace {
public:
  int id = -1;
  int X_ = -1, Y_ = -1;

  GridHalfEdge *lHalfEdge = NULL, *rHalfEdge = NULL, *dHalfEdge = NULL,
               *uHalfEdge = NULL;
  GridVertex *ldVertex = NULL, *rdVertex = NULL, *luVertex = NULL,
             *ruVertex = NULL;

  GridFace *parent = NULL;
  GridFace *ld_son = NULL, *lu_son = NULL, *ru_son = NULL, *rd_son = NULL;

  enum LocateInParent { ld, rd, lu, ru, undefined };
  LocateInParent which_son = undefined;

  enum STATE { in_mesh, outof_mesh, on_boundary, unknown };
  STATE state = outof_mesh;

  bool is_valid = false;
};

class GridHalfEdge {
public:
  int id = -1;
  GridVertex *EndVertex = NULL, *BeginVertex = NULL;
  GridHalfEdge *NextHalfEdge = NULL, *OppositeHalfEdge = NULL;
  GridFace *Face = NULL;
  GridEdge *Edge = NULL;

  GridHalfEdge *Parent = NULL;

  enum E_State { intersect, no_intersect };
  E_State state = no_intersect;
};

class GridEdge {
public:
  bool isHorizontal;
  int id = -1;
  GridHalfEdge *HalfEdge0 = NULL, *HalfEdge1 = NULL;

  enum E_State { intersect, no_intersect };
  E_State state = no_intersect;
};

class GridVertex {
public:
  int id = -1;
  int X = -1, Y = -1;
  std::array<double, 2> coord;

  GridHalfEdge *lBeginHalfEdge = NULL, *lEndHalfEdge = NULL,
               *rBeginHalfEdge = NULL, *rEndHalfEdge = NULL,
               *dBeginHalfEdge = NULL, *dEndHalfEdge = NULL,
               *uBeginHalfEdge = NULL, *uEndHalfEdge = NULL;
  GridFace *ldFace = NULL, *rdFace = NULL, *luFace = NULL, *ruFace = NULL;
  GridEdge *lEdge = NULL, *rEdge = NULL, *uEdge = NULL, *dEdge = NULL;

  GridVertex *parent = NULL;
  GridVertex *son = NULL;

  enum STATE { in_mesh, outof_mesh, on_boundary, unknown };
  STATE state = unknown;
};

class HierarchicalPixelGrid {
public:
  explicit HierarchicalPixelGrid(int grid_size);
  void SetGridBBox(const std::array<double, 2> &bbMin,
                   const std::array<double, 2> &bbMax);
  void Run();
  void PrintFindVEF(std::ofstream &file);
  void PrintQuadMeshOBJ(std::ofstream &file);
  std::vector<GridVertex> V;
  std::vector<GridFace> F;
  std::vector<GridHalfEdge> He;
  std::vector<GridEdge> E;

private:
  void create();
  const int grid_size_;
  std::array<double, 2> bbMin_, bbMax_;
};