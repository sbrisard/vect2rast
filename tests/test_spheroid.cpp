#include <iostream>
#include <numeric>

#include "test_utils.hpp"
#include "vect2rast/spheroid.hpp"

void test_spheroid_belongs(
    const std::array<double, vect2rast::Spheroid::dim>& center, double a,
    double c) {
  constexpr size_t dim = vect2rast::Spheroid::dim;
  std::cout << "test_spheroid_belongs(center="
            << vect2rast::repr(center.cbegin(), center.cend()) << ", a=" << a
            << ", c =" << c << ")...";

  auto directions = generate_directions<dim>();

  double alpha_in = 0.9999;
  double alpha_out = 1.0001;
  std::array<double, 3> ex{1., 0., 0.};
  std::array<double, dim> p_in, p_out;
  for (const auto d3: directions) {
    vect2rast::Spheroid spheroid{center, a, c, d3};
    auto d1 = v2r_cross(ex, d3);
    v2r_normalize(d1);
    auto d2 = v2r_cross(d3, d1);
    for (const auto n: directions) {
      double n1 = std::inner_product(d1.cbegin(), d1.cend(), n.cbegin(), 0.);
      double n2 = std::inner_product(d2.cbegin(), d2.cend(), n.cbegin(), 0.);
      double n3 = std::inner_product(d3.cbegin(), d3.cend(), n.cbegin(), 0.);
      for (size_t k = 0; k < dim; k++) {
        double r = a * (n1 * d1[k] + n2 * d2[k]) + c * n3 * d3[k];
        p_in[k] = spheroid.center[k] + alpha_in * r;
        p_out[k] = spheroid.center[k] + alpha_out * r;
      }
      assert_true(spheroid.belongs(p_in));
      assert_false(spheroid.belongs(p_out));
    }
  }
  std::cout << " OK" << std::endl;
}

void test_spheroid_all() {
  std::array center{1.2, -3.4, 5.6};
  test_spheroid_belongs(center, 0.5, 0.02);
  test_spheroid_belongs(center, 0.02, 0.5);
}
