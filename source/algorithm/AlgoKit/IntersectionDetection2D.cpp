#include "IntersectionDetection2D.h"
namespace AlgoKit {
int Orient2D(const double &px, const double &py,   //
             const double &p0x, const double &p0y, //
             const double &p1x, const double &p1y) {
  // 判定点线关系, 返回-1, 0, 1
  // 0表示共线
  const double dl = (p0x - px) * (p1y - py);
  const double dr = (p0y - py) * (p1x - px);
  const double det = dl - dr;
  const double eb = 3.3306690738754706e-016 * (fabs(dl) + fabs(dr));
  return ((det >= eb) - (-det >= eb));
}

bool CheckPointInTriangle2D(const double &px, const double &py,   //
                            const double &p0x, const double &p0y, //
                            const double &p1x, const double &p1y, //
                            const double &p2x, const double &p2y) {
  double det01 = Orient2D(px, py, p0x, p0y, p1x, p1y);
  double det12 = Orient2D(px, py, p1x, p1y, p2x, p2y);
  double det20 = Orient2D(px, py, p2x, p2y, p0x, p0y);
  return ((det01 >= 0) && (det12 >= 0) && (det20 >= 0)) ||
         ((det01 <= 0) && (det12 <= 0) && (det20 <= 0));
}

bool CheckBoxTriangleIntersection(const std::array<double, 2> &bbmin,
                                  const std::array<double, 2> &bbmax,
                                  const std::array<double, 2> &v0, //
                                  const std::array<double, 2> &v1, //
                                  const std::array<double, 2> &v2) {
  if (v0[0] <= bbmin[0] && v1[0] <= bbmin[0] && v2[0] <= bbmin[0])
    return false;
  if (v0[0] >= bbmax[0] && v1[0] >= bbmax[0] && v2[0] >= bbmax[0])
    return false;
  if (v0[1] <= bbmin[1] && v1[1] <= bbmin[1] && v2[1] <= bbmin[1])
    return false;
  if (v0[1] >= bbmax[1] && v1[1] >= bbmax[1] && v2[1] >= bbmax[1])
    return false;

  const double &x0 = v0[0];
  const double &x1 = v1[0];
  const double &x2 = v2[0];
  const double &y0 = v0[1];
  const double &y1 = v1[1];
  const double &y2 = v2[1];

  const double &xmin = bbmin[0];
  const double &xmax = bbmax[0];
  const double &ymin = bbmin[1];
  const double &ymax = bbmax[1];

  // Consider a triangle on 2d subspace.
  // We use three segments from triangle as separating line.
  // For each separating line, if it separates the third point of triangle and
  // box on 2d subspace, triangle do not intersect with box, as well as in 3d
  // space.

  // So, we test whether 2d triangle intersects with box on x-y, y-z, x-z
  // subspace. If they do not intersect in any subspace, they do not intersect
  // in 3d space. Otherwise, triangle may intersect with box. We need to do
  // triangle-triangle intersection test.

#define Separate2D(seg_x0, seg_y0, seg_x1, seg_y1, px, py, end)                \
  int seg_x_diff_sign = seg_x1 > seg_x0 ? 1 : seg_x1 == seg_x0 ? 0 : -1;       \
  int seg_y_diff_sign = seg_y1 > seg_y0 ? 1 : seg_y1 == seg_y0 ? 0 : -1;       \
  if (seg_x_diff_sign * seg_y_diff_sign == 1) {                                \
    int ori = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, px, py);                \
    if (ori == 0) {                                                            \
      int ori0 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymax);         \
      int ori1 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymin);         \
      if (ori0 * ori1 == 1)                                                    \
        return false;                                                          \
    } else {                                                                   \
      int ori0 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymax);         \
      if (ori0 * ori >= 0)                                                     \
        goto end;                                                              \
      int ori1 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymin);         \
      if (ori1 * ori >= 0)                                                     \
        goto end;                                                              \
      return false;                                                            \
    }                                                                          \
  } else if (seg_x_diff_sign * seg_y_diff_sign == -1) {                        \
    int ori = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, px, py);                \
    if (ori == 0) {                                                            \
      int ori0 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymin);         \
      int ori1 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymax);         \
      if (ori0 * ori1 == 1)                                                    \
        return false;                                                          \
    } else {                                                                   \
      int ori0 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmin, ymin);         \
      if (ori0 * ori >= 0)                                                     \
        goto end;                                                              \
      int ori1 = Orient2D(seg_x0, seg_y0, seg_x1, seg_y1, xmax, ymax);         \
      if (ori1 * ori >= 0)                                                     \
        goto end;                                                              \
      return false;                                                            \
    }                                                                          \
  } else if (seg_x_diff_sign == 0 && seg_y_diff_sign != 0) {                   \
    int ori = px > seg_x0 ? 1 : px == seg_x0 ? 0 : -1;                         \
    if (ori == 0) {                                                            \
      int ori0 = xmin > seg_x0 ? 1 : xmin == seg_x0 ? 0 : -1;                  \
      int ori1 = xmax > seg_x0 ? 1 : xmax == seg_x0 ? 0 : -1;                  \
      if (ori0 * ori1 == 1)                                                    \
        return false;                                                          \
    } else {                                                                   \
      int ori0 = xmin > seg_x0 ? 1 : xmin == seg_x0 ? 0 : -1;                  \
      if (ori0 * ori >= 0)                                                     \
        goto end;                                                              \
      int ori1 = xmax > seg_x0 ? 1 : xmax == seg_x0 ? 0 : -1;                  \
      if (ori1 * ori >= 0)                                                     \
        goto end;                                                              \
      return false;                                                            \
    }                                                                          \
  } else if (seg_y_diff_sign == 0 && seg_x_diff_sign != 0) {                   \
    int ori = py > seg_y0 ? 1 : py == seg_y0 ? 0 : -1;                         \
    if (ori == 0) {                                                            \
      int ori0 = ymin > seg_y0 ? 1 : ymin == seg_y0 ? 0 : -1;                  \
      int ori1 = ymax > seg_y0 ? 1 : ymax == seg_y0 ? 0 : -1;                  \
      if (ori0 * ori1 == 1)                                                    \
        return false;                                                          \
    } else {                                                                   \
      int ori0 = ymin > seg_y0 ? 1 : ymin == seg_y0 ? 0 : -1;                  \
      if (ori0 * ori >= 0)                                                     \
        goto end;                                                              \
      int ori1 = ymax > seg_y0 ? 1 : ymax == seg_y0 ? 0 : -1;                  \
      if (ori1 * ori >= 0)                                                     \
        goto end;                                                              \
      return false;                                                            \
    }                                                                          \
  }

  {
    Separate2D(x0, y0, x1, y1, x2, y2, END_x_y_01);
  }
END_x_y_01:;
  {
    Separate2D(x1, y1, x2, y2, x0, y0, END_x_y_12);
  }
END_x_y_12:;
  {
    Separate2D(x2, y2, x0, y0, x1, y1, END_x_y_20);
  }
END_x_y_20:;

  return true;
}

// 检测box和线段的相交
bool CheckBoxLineSegmentIntersection(const std::array<double, 2> &bbmin,
                                     const std::array<double, 2> &bbmax,
                                     const std::array<double, 2> &v0, //
                                     const std::array<double, 2> &v1) {
  // v0和v1的坐标同时小于或同时大于bbox
  if (v0[0] < bbmin[0] && v1[0] < bbmin[0])
    return false;
  if (v0[0] > bbmax[0] && v1[0] > bbmax[0])
    return false;
  if (v0[1] < bbmin[1] && v1[1] < bbmin[1])
    return false;
  if (v0[1] > bbmax[1] && v1[1] > bbmax[1])
    return false;

  // v0或v1位于bbox内部
  if (v0[0] >= bbmin[0] && v0[0] <= bbmax[0] && v0[1] >= bbmin[1] &&
      v0[1] <= bbmax[1])
    return true;
  if (v1[0] >= bbmin[0] && v1[0] <= bbmax[0] && v1[1] >= bbmin[1] &&
      v1[1] <= bbmax[1])
    return true;

  int x_diff_sign = v1[0] > v0[0] ? 1 : v1[0] == v0[0] ? 0 : -1;
  int y_diff_sign = v1[1] > v0[1] ? 1 : v1[1] == v0[1] ? 0 : -1;

  if (x_diff_sign * y_diff_sign == 1) {
    // seg向上倾斜, 算bbox的lu点和rd点与seg的关系
    int ori_bbox_left_up =
        Orient2D(bbmin[0], bbmax[1], v0[0], v0[1], v1[0], v1[1]);
    if (ori_bbox_left_up == 0)
      return true;
    int ori_bbox_right_down =
        Orient2D(bbmax[0], bbmin[1], v0[0], v0[1], v1[0], v1[1]);
    if (ori_bbox_right_down == 0)
      return true;
    return ori_bbox_right_down * ori_bbox_left_up == -1;
  } else if (x_diff_sign * y_diff_sign == -1) {
    // seg向下倾斜, 算bbox的ld点和ru点与seg的关系
    int ori_bbox_left_down =
        Orient2D(bbmin[0], bbmin[1], v0[0], v0[1], v1[0], v1[1]);
    if (ori_bbox_left_down == 0)
      return true;
    int ori_bbox_right_up =
        Orient2D(bbmax[0], bbmax[1], v0[0], v0[1], v1[0], v1[1]);
    if (ori_bbox_right_up == 0)
      return true;
    return ori_bbox_right_up * ori_bbox_left_down == -1;
  } else if (x_diff_sign == 0 && y_diff_sign != 0) {
    // seg竖直, 算bbox的lu点和rd点与seg的关系
    return v0[0] >= bbmin[0] && v0[0] <= bbmax[0];
  } else if (x_diff_sign != 0 && y_diff_sign == 0) {
    // seg水平, 算bbox的ld点和ru点与seg的关系
    return v0[1] >= bbmin[1] && v0[1] <= bbmax[1];
  } else if (x_diff_sign == 0 && y_diff_sign == 0)
    // seg退化为单个点, 由于前面检测过v0/v1是否在内部, 此处可以直接判负
    return false;
  return true;
}
} // namespace AlgoKit