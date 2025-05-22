#include "HierarchicalPixelGrid.h"

HierarchicalPixelGrid::HierarchicalPixelGrid(int grid_size)
    : grid_size_(grid_size) {}

void HierarchicalPixelGrid::create() {
  double pixel_unit_size = 1.0 / grid_size_;

  V.resize((grid_size_ + 1) * (grid_size_ + 1));
  for (int y_ = 0; y_ < grid_size_ + 1; y_++) {
    for (int x_ = 0; x_ < grid_size_ + 1; x_++) {
      GridVertex &v = V[y_ * (grid_size_ + 1) + x_];
      v.X = x_;
      v.Y = y_;
      v.id = x_ + y_ * (grid_size_ + 1);
      v.coord = {x_ * pixel_unit_size, y_ * pixel_unit_size};
    }
  }

  F.resize(grid_size_ * grid_size_);
  for (int y_ = 0; y_ < grid_size_; ++y_) {
    for (int x_ = 0; x_ < grid_size_; ++x_) {
      GridFace &f = F[y_ * grid_size_ + x_];
      f.X_ = x_;
      f.Y_ = y_;
      f.id = x_ + y_ * grid_size_;
      f.ldVertex = &V[y_ * (grid_size_ + 1) + x_];
      V[y_ * (grid_size_ + 1) + x_].ruFace = &f;

      f.rdVertex = &V[y_ * (grid_size_ + 1) + x_ + 1];
      V[y_ * (grid_size_ + 1) + x_ + 1].luFace = &f;

      f.luVertex = &V[(y_ + 1) * (grid_size_ + 1) + x_];
      V[(y_ + 1) * (grid_size_ + 1) + x_].rdFace = &f;

      f.ruVertex = &V[(y_ + 1) * (grid_size_ + 1) + x_ + 1];
      V[(y_ + 1) * (grid_size_ + 1) + x_ + 1].ldFace = &f;
    }
  }

  He.resize(4 * grid_size_ * grid_size_ + 2 * (grid_size_ + grid_size_));

  for (int i = 0; i < static_cast<int>(F.size()); ++i) {
    GridHalfEdge &he0 = He[4 * i + 0];
    GridHalfEdge &he1 = He[4 * i + 1];
    GridHalfEdge &he2 = He[4 * i + 2];
    GridHalfEdge &he3 = He[4 * i + 3];

    he0.id = 4 * i + 0;
    he1.id = 4 * i + 1;
    he2.id = 4 * i + 2;
    he3.id = 4 * i + 3;

    he0.BeginVertex = F[i].rdVertex;
    F[i].rdVertex->lBeginHalfEdge = &he0;
    he0.EndVertex = F[i].ldVertex;
    F[i].ldVertex->rEndHalfEdge = &he0;
    he0.Face = &F[i];
    F[i].dHalfEdge = &he0;
    he0.NextHalfEdge = &he1;

    he1.BeginVertex = F[i].ldVertex;
    F[i].ldVertex->uBeginHalfEdge = &he1;
    he1.EndVertex = F[i].luVertex;
    F[i].luVertex->dEndHalfEdge = &he1;
    he1.Face = &F[i];
    F[i].lHalfEdge = &he1;
    he1.NextHalfEdge = &he2;

    he2.BeginVertex = F[i].luVertex;
    F[i].luVertex->rBeginHalfEdge = &he2;
    he2.EndVertex = F[i].ruVertex;
    F[i].ruVertex->lEndHalfEdge = &he2;
    he2.Face = &F[i];
    F[i].uHalfEdge = &he2;
    he2.NextHalfEdge = &he3;

    he3.BeginVertex = F[i].ruVertex;
    F[i].ruVertex->dBeginHalfEdge = &he3;
    he3.EndVertex = F[i].rdVertex;
    F[i].rdVertex->uEndHalfEdge = &he3;
    he3.Face = &F[i];
    F[i].rHalfEdge = &he3;
    he3.NextHalfEdge = &he0;
  }

  E.resize(2 * grid_size_ * grid_size_ + grid_size_ + grid_size_);
  for (int y_ = 0; y_ < grid_size_ + 1; ++y_)
    for (int x_ = 0; x_ < grid_size_; ++x_) {
      GridEdge &e = E[y_ * grid_size_ + x_];
      e.id = y_ * grid_size_ + x_;
      e.isHorizontal = true;

      GridVertex &v0 = V[y_ * (grid_size_ + 1) + x_];
      GridVertex &v1 = V[y_ * (grid_size_ + 1) + x_ + 1];

      v0.rEdge = &e;
      v1.lEdge = &e;

      if (v0.rBeginHalfEdge != NULL && v1.lBeginHalfEdge != NULL) {
        e.HalfEdge0 = v0.rBeginHalfEdge;
        e.HalfEdge1 = v1.lBeginHalfEdge;

        v0.rBeginHalfEdge->OppositeHalfEdge = v1.lBeginHalfEdge;
        v0.rBeginHalfEdge->Edge = &e;

        v1.lBeginHalfEdge->OppositeHalfEdge = v0.rBeginHalfEdge;
        v1.lBeginHalfEdge->Edge = &e;
      } else if (v0.rBeginHalfEdge != NULL && v1.lBeginHalfEdge == NULL) {
        e.HalfEdge0 = v0.rBeginHalfEdge;
        v0.rBeginHalfEdge->Edge = &e;
      } else if (v0.rBeginHalfEdge == NULL && v1.lBeginHalfEdge != NULL) {
        e.HalfEdge1 = v1.lBeginHalfEdge;
        v1.lBeginHalfEdge->Edge = &e;
      }
    }

  int baseE = grid_size_ * (grid_size_ + 1);
  for (int y_ = 0; y_ < grid_size_; ++y_)
    for (int x_ = 0; x_ < grid_size_ + 1; ++x_) {
      GridEdge &e = E[baseE + y_ * (grid_size_ + 1) + x_];
      e.id = baseE + y_ * (grid_size_ + 1) + x_;
      e.isHorizontal = false;

      GridVertex &v0 = V[y_ * (grid_size_ + 1) + x_];
      GridVertex &v1 = V[(y_ + 1) * (grid_size_ + 1) + x_];

      v0.uEdge = &e;
      v1.dEdge = &e;

      if (v0.uBeginHalfEdge != NULL && v1.dBeginHalfEdge != NULL) {
        e.HalfEdge0 = v0.uBeginHalfEdge;
        e.HalfEdge1 = v1.dBeginHalfEdge;

        v0.uBeginHalfEdge->OppositeHalfEdge = v1.dBeginHalfEdge;
        v0.uBeginHalfEdge->Edge = &e;

        v1.dBeginHalfEdge->OppositeHalfEdge = v0.uBeginHalfEdge;
        v1.dBeginHalfEdge->Edge = &e;
      } else if (v0.uBeginHalfEdge != NULL && v1.dBeginHalfEdge == NULL) {
        e.HalfEdge0 = v0.uBeginHalfEdge;
        v0.uBeginHalfEdge->Edge = &e;
      } else if (v0.uBeginHalfEdge == NULL && v1.dBeginHalfEdge != NULL) {
        e.HalfEdge1 = v1.dBeginHalfEdge;
        v1.dBeginHalfEdge->Edge = &e;
      }
    }

  int HeId = 4 * grid_size_ * grid_size_;
  for (int i = 0; i < 4 * grid_size_ * grid_size_; i++) {
    GridHalfEdge &he0 = He[i];
    if (he0.OppositeHalfEdge != NULL)
      continue;
    GridHalfEdge &he1 = He[HeId];
    he1.id = HeId;
    // he1.BeginVertex = he0.EndVertex;
    // he1.EndVertex = he0.BeginVertex;
    he1.OppositeHalfEdge = &he0;
    he0.OppositeHalfEdge = &he1;
    he1.Edge = he0.Edge;
    if (he0.Edge->HalfEdge0 == NULL)
      he0.Edge->HalfEdge0 = &he1;
    else
      he0.Edge->HalfEdge1 = &he1;
    HeId++;
  }
}
