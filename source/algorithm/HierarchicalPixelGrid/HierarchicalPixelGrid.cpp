#include "HierarchicalPixelGrid.h"

HierarchicalPixelGrid::HierarchicalPixelGrid(int grid_size)
    : grid_size_(grid_size) {
  create();
}

void HierarchicalPixelGrid::create() {

  V.resize((grid_size_ + 1) * (grid_size_ + 1));
  for (int y_ = 0; y_ < grid_size_ + 1; y_++) {
    for (int x_ = 0; x_ < grid_size_ + 1; x_++) {
      GridVertex &v = V[y_ * (grid_size_ + 1) + x_];
      v.X = x_;
      v.Y = y_;
      v.id = x_ + y_ * (grid_size_ + 1);
      v.coord = {x_, y_};
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

void HierarchicalPixelGrid::PrintFindVEF(std::ofstream &file) {
  file << "PE" << std::endl;
  for (std::size_t i = 0; i < F.size(); ++i) {
    if (F[i].state != GridFace::in_mesh)
      continue;
    GridVertex *vertex_0 = F[i].ldVertex;
    GridVertex *vertex_1 = F[i].luVertex;
    GridVertex *vertex_2 = F[i].ruVertex;
    GridVertex *vertex_3 = F[i].rdVertex;
    file << std::format("{} {} 0 {} {} 0 0 0 64", //
                        vertex_0->coord[0] / grid_size_,
                        vertex_0->coord[1] / grid_size_,
                        vertex_1->coord[0] / grid_size_,
                        vertex_1->coord[1] / grid_size_)
         << std::endl;
    file << std::format("{} {} 0 {} {} 0 0 0 64", //
                        vertex_1->coord[0] / grid_size_,
                        vertex_1->coord[1] / grid_size_,
                        vertex_2->coord[0] / grid_size_,
                        vertex_2->coord[1] / grid_size_)
         << std::endl;
    file << std::format("{} {} 0 {} {} 0 0 0 64", //
                        vertex_2->coord[0] / grid_size_,
                        vertex_2->coord[1] / grid_size_,
                        vertex_3->coord[0] / grid_size_,
                        vertex_3->coord[1] / grid_size_)
         << std::endl;
    file << std::format("{} {} 0 {} {} 0 0 0 64", //
                        vertex_3->coord[0] / grid_size_,
                        vertex_3->coord[1] / grid_size_,
                        vertex_0->coord[0] / grid_size_,
                        vertex_0->coord[1] / grid_size_)
         << std::endl;
  }
}

void HierarchicalPixelGrid::PrintQuadMeshOBJ(std::ofstream &file) {
  std::vector<int> map_grid_vid_to_quad_vid(V.size(), -1);
  std::vector<int> map_quad_vid_to_grid_vid;
  map_quad_vid_to_grid_vid.reserve(V.size());
  std::size_t quad_idx_ = 0;
  for (std::size_t i = 0; i < F.size(); ++i) {
    if (F[i].state != GridFace::in_mesh)
      continue;
    GridVertex *vertex_0 = F[i].ldVertex;
    GridVertex *vertex_1 = F[i].luVertex;
    GridVertex *vertex_2 = F[i].ruVertex;
    GridVertex *vertex_3 = F[i].rdVertex;

    if (map_grid_vid_to_quad_vid[vertex_0->id] == -1) {
      map_grid_vid_to_quad_vid[static_cast<std::size_t>(vertex_0->id)] =
          static_cast<int>(quad_idx_);
      map_quad_vid_to_grid_vid.push_back(vertex_0->id);
      quad_idx_++;
    }
    if (map_grid_vid_to_quad_vid[vertex_1->id] == -1) {
      map_grid_vid_to_quad_vid[static_cast<std::size_t>(vertex_1->id)] =
          static_cast<int>(quad_idx_);
      map_quad_vid_to_grid_vid.push_back(vertex_1->id);
      quad_idx_++;
    }
    if (map_grid_vid_to_quad_vid[vertex_2->id] == -1) {
      map_grid_vid_to_quad_vid[static_cast<std::size_t>(vertex_2->id)] =
          static_cast<int>(quad_idx_);
      map_quad_vid_to_grid_vid.push_back(vertex_2->id);
      quad_idx_++;
    }
    if (map_grid_vid_to_quad_vid[vertex_3->id] == -1) {
      map_grid_vid_to_quad_vid[static_cast<std::size_t>(vertex_3->id)] =
          static_cast<int>(quad_idx_);
      map_quad_vid_to_grid_vid.push_back(vertex_3->id);
      quad_idx_++;
    }
  }

  for (std::size_t i = 0; i < quad_idx_; ++i) {
    file << std::format("v {} {} 0",
                        V[map_quad_vid_to_grid_vid[i]].coord[0] / grid_size_,
                        V[map_quad_vid_to_grid_vid[i]].coord[1] / grid_size_)
         << std::endl;
  }

  for (std::size_t i = 0; i < F.size(); ++i) {
    if (F[i].state != GridFace::in_mesh)
      continue;
    file << std::format("f {} {} {} {}",
                        map_grid_vid_to_quad_vid[F[i].ldVertex->id] + 1,
                        map_grid_vid_to_quad_vid[F[i].rdVertex->id] + 1,
                        map_grid_vid_to_quad_vid[F[i].ruVertex->id] + 1,
                        map_grid_vid_to_quad_vid[F[i].luVertex->id] + 1)
         << std::endl;
  }
}

GridElement
HierarchicalPixelGrid::LocateGridElement(const std::array<double, 2> &coord) {
  constexpr double EPS = 1e-6;

  if (coord[0] < 0.0 - EPS || coord[0] > 1.0 + EPS || coord[1] < 0.0 - EPS ||
      coord[1] > 1.0 + EPS)
    return std::monostate();

  const double x_scaled = coord[0] * grid_size_;
  const double y_scaled = coord[1] * grid_size_;

  // 3. 分离整数和小数部分 (考虑边界)
  auto decompose = [&](double val) -> std::pair<int, double> {
    int integer = static_cast<int>(std::floor(val + EPS)); // 防止0.999999被误判
    double fraction = val - integer;
    // 处理右边界 (e.g. grid_size=5, val=5.0 → integer=5)
    if (std::abs(fraction - 1.0) < EPS) {
      integer += 1;
      fraction = 0.0;
    }
    return {integer, fraction};
  };
  const std::pair<int, double> x_decompose = decompose(x_scaled);
  const std::pair<int, double> y_decompose = decompose(y_scaled);

  // 4. 判断是否为顶点 (双坐标均为整数)
  const bool x_is_int =
      (x_decompose.second < EPS) || (1.0 - x_decompose.second < EPS);
  const bool y_is_int =
      (y_decompose.second < EPS) || (1.0 - y_decompose.second < EPS);

  // 5. 判断是否为边 (x坐标为整数，y坐标为小数)
  if (x_is_int && y_is_int) {
    return V[y_decompose.first * (grid_size_ + 1) + x_decompose.first];
  } else if (!x_is_int && y_is_int) {
    GridVertex &vertex =
        V[y_decompose.first * (grid_size_ + 1) + x_decompose.first];
    return (vertex.rBeginHalfEdge != NULL) ? (*vertex.rBeginHalfEdge->Edge)
                                           : (*vertex.rEndHalfEdge->Edge);
  } else if (x_is_int && !y_is_int) {
    GridVertex &vertex =
        V[y_decompose.first * (grid_size_ + 1) + x_decompose.first];
    return (vertex.uBeginHalfEdge != NULL) ? (*vertex.uBeginHalfEdge->Edge)
                                           : (*vertex.uEndHalfEdge->Edge);
  } else {
    return F[y_decompose.first * grid_size_ + x_decompose.first];
  }

  return std::monostate();
}