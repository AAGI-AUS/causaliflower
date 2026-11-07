// =============================================================================
// native.cpp


// Native projection and curve contact steps with complete R fallbacks
// =============================================================================


// These functions implement the package's R algorithms in C++.
// No dagitty code is incorporated. R retains vertex sorting, edge ordering,
// curve coordinates, candidate construction and final paths.
#define R_NO_REMAP
#include <R.h>
#include <Rinternals.h>
#include <R_ext/Utils.h>
#include <R_ext/Rdynload.h>
#include <R_ext/Visibility.h>
#include <algorithm>
#include <cmath>
#include <cstdint>
#include <cstring>
#include <limits>


extern "C" SEXP cf_project_edges(SEXP count, SEXP from, SEXP to, SEXP arc, SEXP omitted) { // Project edges through omitted vertices
  if (TYPEOF(count) != INTSXP || XLENGTH(count) != 1 ||
      INTEGER(count)[0] < 0 || TYPEOF(from) != INTSXP ||
      TYPEOF(to) != INTSXP || TYPEOF(arc) != LGLSXP ||
      TYPEOF(omitted) != LGLSXP) {
    Rf_error("Invalid native projection inputs.");
  }
  const int node_count = INTEGER(count)[0];
  const R_xlen_t edge_count = XLENGTH(from);
  if (XLENGTH(to) != edge_count || XLENGTH(arc) != edge_count ||
      XLENGTH(omitted) != node_count) {
    Rf_error("Native projection input lengths differ.");
  }
  if (static_cast<size_t>(node_count) > std::numeric_limits<size_t>::max() /
      (node_count ? static_cast<size_t>(node_count) : 1)) {
    Rf_error("Native projection matrix is too large.");
  }
  const size_t matrix_size = static_cast<size_t>(node_count) * node_count;
  // Temporary storage allocated by R is reclaimed even if an interrupt unwinds this
  // call, unlike C++ containers skipped by R's nonlocal return.
  auto *directed = reinterpret_cast<uint8_t *>(R_alloc(matrix_size, sizeof(uint8_t)));
  auto *bidirected = reinterpret_cast<uint8_t *>(R_alloc(matrix_size, sizeof(uint8_t)));
  std::memset(directed, 0, matrix_size);
  std::memset(bidirected, 0, matrix_size);
  const auto position = [node_count](int source, int target) { // Index adjacency arrays stored by row
    return static_cast<size_t>(source) * node_count + target;
  };
  for (R_xlen_t edge = 0; edge < edge_count; ++edge) {
    const int source = INTEGER(from)[edge] - 1;
    const int target = INTEGER(to)[edge] - 1;
    if (source < 0 || source >= node_count || target < 0 || target >= node_count ||
        LOGICAL(arc)[edge] == NA_LOGICAL) {
      Rf_error("Invalid native projection edge.");
    }
    if (LOGICAL(arc)[edge]) {
      bidirected[position(source, target)] = 1;
      bidirected[position(target, source)] = 1;
    } else {
      directed[position(source, target)] = 1;
    }
  }
  auto *parents = reinterpret_cast<int *>(R_alloc(node_count, sizeof(int)));
  auto *children = reinterpret_cast<int *>(R_alloc(node_count, sizeof(int)));
  auto *spouses = reinterpret_cast<int *>(R_alloc(node_count, sizeof(int)));
  for (int node = 0; node < node_count; ++node) {
    if (LOGICAL(omitted)[node] == NA_LOGICAL) {
      Rf_error("Invalid native projection omission.");
    }
    if (!LOGICAL(omitted)[node]) {
      continue;
    }
    R_CheckUserInterrupt();
    int parent_count = 0;
    int child_count = 0;
    int spouse_count = 0;
    for (int other = 0; other < node_count; ++other) {
      if (directed[position(other, node)]) {
        parents[parent_count++] = other;
      }
      if (directed[position(node, other)]) {
        children[child_count++] = other;
      }
      if (bidirected[position(node, other)]) {
        spouses[spouse_count++] = other;
      }
    }
    for (int parent = 0; parent < parent_count; ++parent) {
      for (int child = 0; child < child_count; ++child) {
        directed[position(parents[parent], children[child])] = 1;
      }
    }
    if (child_count > 1) {
      for (int left = 0; left < child_count; ++left) {
        for (int right = 0; right < child_count; ++right) {
          bidirected[position(children[left], children[right])] = 1;
        }
      }
    }
    for (int spouse = 0; spouse < spouse_count; ++spouse) {
      for (int child = 0; child < child_count; ++child) {
        bidirected[position(spouses[spouse], children[child])] = 1;
        bidirected[position(children[child], spouses[spouse])] = 1;
      }
    }
    for (int other = 0; other < node_count; ++other) {
      bidirected[position(other, other)] = 0;
      directed[position(node, other)] = directed[position(other, node)] = 0;
      bidirected[position(node, other)] = bidirected[position(other, node)] = 0;
    }
  }
  size_t rows = 0;
  for (int target = 0; target < node_count; ++target) {
    for (int source = 0; source < node_count; ++source) {
      rows += directed[position(source, target)] +
        (source < target && bidirected[position(source, target)]);
    }
  }
  if (rows > static_cast<size_t>(std::numeric_limits<int>::max())) {
    Rf_error("Native projection result is too large.");
  }
  SEXP result = PROTECT(Rf_allocMatrix(INTSXP, static_cast<int>(rows), 3));
  size_t row = 0;
  for (int type = 0; type < 2; ++type) {
    for (int target = 0; target < node_count; ++target) {
      for (int source = 0; source < node_count; ++source) {
        if ((!type && directed[position(source, target)]) ||
            (type && source < target && bidirected[position(source, target)])) {
          INTEGER(result)[row] = source + 1;
          INTEGER(result)[row + rows] = target + 1;
          INTEGER(result)[row + 2 * rows] = type;
          ++row;
        }
      }
    }
  }
  UNPROTECT(1);
  return result;
}


namespace {

SEXP field(SEXP value, const char *name) { // Read one named list field without changing its R representation
  if (TYPEOF(value) != VECSXP) {
    return R_NilValue;
  }
  SEXP names = Rf_getAttrib(value, R_NamesSymbol);
  for (R_xlen_t index = 0; index < XLENGTH(value); ++index) {
    if (names != R_NilValue && std::strcmp(CHAR(STRING_ELT(names, index)), name) == 0) {
      return VECTOR_ELT(value, index);
    }
  }
  return R_NilValue;
}


bool number_vector(SEXP value, R_xlen_t length) { // Use the R function for extreme or nonfinite inputs
  if (TYPEOF(value) != REALSXP || XLENGTH(value) != length) {
    return false;
  }
  const double maximum_coordinate = 1e100; // Ordinary planner coordinates are normalised around [0, 1]
  for (R_xlen_t index = 0; index < length; ++index) {
    if (!R_FINITE(REAL(value)[index]) || std::abs(REAL(value)[index]) > maximum_coordinate) {
      return false;
    }
  }
  return true;
}


const char *string(SEXP value) { // Read one UTF-8 string that is not missing
  if (TYPEOF(value) != STRSXP || XLENGTH(value) != 1 || STRING_ELT(value, 0) == NA_STRING) {
    return nullptr;
  }
  return Rf_translateCharUTF8(STRING_ELT(value, 0));
}


struct CurveState {
  const double *start = nullptr, *end = nullptr, *low = nullptr, *high = nullptr;
  const double *points = nullptr, *segment_low = nullptr, *segment_high = nullptr;
  const char *from = nullptr, *to = nullptr;
  double height = 0, sign = 0;
  int rows = 0;
  bool lateral = false;
};


bool read_state(SEXP value, CurveState &state) { // Validate the R planner's cached geometry before using its arrays
  SEXP low = field(value, "low");
  SEXP high = field(value, "high");
  SEXP lateral = field(value, "lateral");
  if (!number_vector(low, 2) || !number_vector(high, 2) ||
      TYPEOF(lateral) != LGLSXP || XLENGTH(lateral) != 1 ||
      LOGICAL(lateral)[0] == NA_LOGICAL) {
    return false;
  }
  state.low = REAL(low);
  state.high = REAL(high);
  state.lateral = LOGICAL(lateral)[0];
  state.from = string(field(value, "from"));
  state.to = string(field(value, "to"));
  if (!state.from || !state.to) {
    return false;
  }
  if (!state.lateral) {
    SEXP start = field(value, "a");
    SEXP end = field(value, "b");
    SEXP height = field(value, "h");
    SEXP sign = field(value, "sign");
    if (!number_vector(start, 2) || !number_vector(end, 2) || !number_vector(height, 1) ||
        !((TYPEOF(sign) == REALSXP && number_vector(sign, 1)) ||
          (TYPEOF(sign) == INTSXP && XLENGTH(sign) == 1 && INTEGER(sign)[0] != NA_INTEGER))) {
      return false;
    }
    state.start = REAL(start);
    state.end = REAL(end);
    state.height = REAL(height)[0];
    state.sign = Rf_asReal(sign);
    const double minimum_span = 1e-100; // Use the R function when a curve that is not lateral is nearly vertical
    if (std::abs(state.end[0] - state.start[0]) < minimum_span) {
      return false;
    }
  }
  SEXP points = field(value, "points");
  SEXP dimensions = Rf_getAttrib(points, R_DimSymbol);
  if (TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2 || INTEGER(dimensions)[1] != 2 ||
      INTEGER(dimensions)[0] < 2 || !number_vector(points, 2 * INTEGER(dimensions)[0])) {
    return false;
  }
  state.rows = INTEGER(dimensions)[0];
  state.points = REAL(points);
  SEXP lower = field(value, "segment_low");
  SEXP upper = field(value, "segment_high");
  if (!number_vector(lower, 2 * (state.rows - 1)) || !number_vector(upper, 2 * (state.rows - 1))) {
    return false;
  }
  state.segment_low = REAL(lower);
  state.segment_high = REAL(upper);
  return true;
}


struct SharedEndpoints {
  double x[2], y[2];
  int size = 0;

  bool contains_point(double point_x, double point_y) const { // Exclude a shared endpoint from a curve contact
    const double squared_tolerance = 1e-14; // Match the R predicate's threshold for squared distance
    for (int index = 0; index < size; ++index) {
      const double delta_x = point_x - x[index];
      const double delta_y = point_y - y[index];
      long double distance = 0; // R sum() accumulates rounded products in long double
      distance += delta_x * delta_x;
      distance += delta_y * delta_y;
      if (static_cast<double>(distance) < squared_tolerance) {
        return true;
      }
    }
    return false;
  }
};


double cross(double first_x, double first_y, double second_x, double second_y) { // Compute a cross product in two dimensions
  return first_x * second_y - first_y * second_x;
}


double sum_pair(double first, double second) { // Match R sum() without changing when the inputs are rounded
  long double sum = 0;
  sum += first;
  sum += second;
  return static_cast<double>(sum);
}


double y_at(const CurveState &state, double x) { // Evaluate the planner's quadratic curve in its original operation order
  const double position = (x - state.start[0]) / (state.end[0] - state.start[0]);
  return state.start[1] + position * (state.end[1] - state.start[1]) +
    state.sign * 4 * state.height * position * (1 - position);
}


bool poly_contact(const CurveState &first, const CurveState &second, const SharedEndpoints &shared) { // Check contacts between sampled path segments
  const double *first_points = first.points;
  const double *second_points = second.points;
  for (int first_row = 0; first_row < first.rows - 1; ++first_row) {
    const double first_x = first_points[first_row];
    const double first_y = first_points[first_row + first.rows];
    const double first_dx = first_points[first_row + 1] - first_x;
    const double first_dy = first_points[first_row + 1 + first.rows] - first_y;
    const double low_x = std::min(first_points[first_row], first_points[first_row + 1]);
    const double low_y = std::min(first_points[first_row + first.rows], first_points[first_row + 1 + first.rows]);
    const double high_x = std::max(first_points[first_row], first_points[first_row + 1]);
    const double high_y = std::max(first_points[first_row + first.rows], first_points[first_row + 1 + first.rows]);
    for (int second_row = 0; second_row < second.rows - 1; ++second_row) {
      if (!(second.segment_low[second_row] <= high_x + 1e-10 &&
            second.segment_high[second_row] >= low_x - 1e-10 &&
            second.segment_low[second_row + second.rows - 1] <= high_y + 1e-10 &&
            second.segment_high[second_row + second.rows - 1] >= low_y - 1e-10)) {
        continue;
      }
      const double second_x = second_points[second_row];
      const double second_y = second_points[second_row + second.rows];
      const double second_dx = second_points[second_row + 1] - second_x;
      const double second_dy = second_points[second_row + 1 + second.rows] - second_y;
      const double denominator = cross(first_dx, first_dy, second_dx, second_dy);
      if (std::abs(denominator) > 1e-12) {
        const double first_position = cross(second_x - first_x, second_y - first_y, second_dx, second_dy) / denominator;
        const double second_position = cross(second_x - first_x, second_y - first_y, first_dx, first_dy) / denominator;
        if (first_position >= -1e-9 && first_position <= 1 + 1e-9 &&
            second_position >= -1e-9 && second_position <= 1 + 1e-9 &&
            !shared.contains_point(first_x + first_position * first_dx, first_y + first_position * first_dy)) {
          return true;
        }
      } else if (sum_pair(first_dx * first_dx, first_dy * first_dy) > 1e-24 &&
                 std::abs(cross(second_x - first_x, second_y - first_y, first_dx, first_dy)) < 1e-10) {
        const double squared_length = sum_pair(first_dx * first_dx, first_dy * first_dy);
        const double start = sum_pair((second_x - first_x) * first_dx, (second_y - first_y) * first_dy) / squared_length;
        const double end = sum_pair((second_x + second_dx - first_x) * first_dx,
                                    (second_y + second_dy - first_y) * first_dy) / squared_length;
        const double low = std::max(0.0, std::min(start, end));
        const double high = std::min(1.0, std::max(start, end));
        if (high - low > 1e-8 ||
            (high >= low && !shared.contains_point(first_x + low * first_dx, first_y + low * first_dy))) {
          return true;
        }
      }
    }
  }
  return false;
}


bool contact(const CurveState &first, const CurveState &second, const SharedEndpoints &shared) { // Check curve intersections with the R predicate's tolerances
  for (int axis = 0; axis < 2; ++axis) {
    if (first.low[axis] > second.high[axis] + 1e-9 || second.low[axis] > first.high[axis] + 1e-9) {
      return false;
    }
  }
  if (first.lateral || second.lateral) {
    return poly_contact(first, second, shared);
  }
  const double low = std::max(std::min(first.start[0], first.end[0]), std::min(second.start[0], second.end[0]));
  const double high = std::min(std::max(first.start[0], first.end[0]), std::max(second.start[0], second.end[0]));
  if (high < low - 1e-10) {
    return false;
  }
  if (high - low < 1e-10) {
    return std::abs(y_at(first, low) - y_at(second, low)) < 1e-9 && !shared.contains_point(low, y_at(first, low));
  }
  const double middle = (low + high) / 2;
  const double difference_left = y_at(first, low) - y_at(second, low);
  const double difference_middle = y_at(first, middle) - y_at(second, middle);
  const double difference_right = y_at(first, high) - y_at(second, high);
  const double quadratic = 2 * (difference_right + difference_left - 2 * difference_middle);
  const double linear = 4 * difference_middle - difference_right - 3 * difference_left;
  const double constant = difference_left;
  const double scale = std::max(std::max(std::abs(quadratic), std::abs(linear)), std::abs(constant));
  if (scale < 1e-10) {
    return true;
  }
  double roots[2];
  int root_count = 0;
  if (std::abs(quadratic) < 1e-10) {
    if (std::abs(linear) >= 1e-10) {
      roots[root_count++] = -constant / linear;
    }
  } else {
    const double discriminant = linear * linear - 4 * quadratic * constant;
    const double tolerance = 1e-12 * (linear * linear + std::abs(4 * quadratic * constant));
    if (discriminant >= -tolerance) {
      const double square_root = std::sqrt(std::max(0.0, discriminant));
      roots[root_count++] = (-linear + -1 * square_root) / (2 * quadratic);
      roots[root_count++] = (-linear + 1 * square_root) / (2 * quadratic);
    }
  }
  for (int index = 0; index < root_count; ++index) {
    const double position = roots[index];
    if (!(position >= -1e-8 && position <= 1 + 1e-8)) {
      continue;
    }
    if (!(std::abs(quadratic * (position * position) + linear * position + constant) <=
          1e-9 * std::max(scale, 1e-12))) {
      continue;
    }
    const double x = low + position * (high - low);
    if (!shared.contains_point(x, y_at(first, x))) {
      return true;
    }
  }
  return false;
}

}


extern "C" SEXP cf_curve_contact(SEXP first, SEXP second, SEXP xy) { // Return NA when the R predicate must handle the supplied geometry
  CurveState first_state, second_state;
  if (!read_state(first, first_state) || !read_state(second, second_state)) {
    return Rf_ScalarLogical(NA_LOGICAL);
  }
  SEXP dimensions = Rf_getAttrib(xy, R_DimSymbol);
  SEXP dimnames = Rf_getAttrib(xy, R_DimNamesSymbol);
  if (TYPEOF(dimensions) != INTSXP || XLENGTH(dimensions) != 2 ||
      INTEGER(dimensions)[1] != 2 || TYPEOF(dimnames) != VECSXP || XLENGTH(dimnames) != 2) {
    return Rf_ScalarLogical(NA_LOGICAL);
  }
  const int node_count = INTEGER(dimensions)[0];
  SEXP names = VECTOR_ELT(dimnames, 0);
  if (TYPEOF(names) != STRSXP || XLENGTH(names) != node_count || !number_vector(xy, 2 * node_count)) {
    return Rf_ScalarLogical(NA_LOGICAL);
  }
  SharedEndpoints shared;
  const char *endpoint_ids[2] = {first_state.from, first_state.to};
  for (int endpoint = 0; endpoint < 2; ++endpoint) {
    if (endpoint && std::strcmp(endpoint_ids[0], endpoint_ids[1]) == 0) {
      continue;
    }
    if (std::strcmp(endpoint_ids[endpoint], second_state.from) &&
        std::strcmp(endpoint_ids[endpoint], second_state.to)) {
      continue;
    }
    bool found = false;
    for (int index = 0; index < node_count; ++index) {
      if (STRING_ELT(names, index) != NA_STRING &&
          !std::strcmp(endpoint_ids[endpoint], Rf_translateCharUTF8(STRING_ELT(names, index)))) {
        shared.x[shared.size] = REAL(xy)[index];
        shared.y[shared.size++] = REAL(xy)[index + node_count];
        found = true;
        break;
      }
    }
    if (!found) {
      return Rf_ScalarLogical(NA_LOGICAL);
    }
  }
  return Rf_ScalarLogical(contact(first_state, second_state, shared));
}


static const R_CallMethodDef call_methods[] = { // Registered entry points and exact argument counts
  {"cf_project_edges", reinterpret_cast<DL_FUNC>(&cf_project_edges), 5},
  {"cf_curve_contact", reinterpret_cast<DL_FUNC>(&cf_curve_contact), 3},
  {nullptr, nullptr, 0}
};


extern "C" void attribute_visible R_init_causaliflower(DllInfo *dll) { // Require registered symbols for both kernels
  R_registerRoutines(dll, nullptr, call_methods, nullptr, nullptr);
  R_useDynamicSymbols(dll, FALSE);
  R_forceSymbols(dll, TRUE);
  return;
}
