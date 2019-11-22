#include <glib.h>

void v2r_setup_test_ndsphere();
void v2r_setup_test_spheroid();

int main(int argc, char **argv) {
  g_test_init(&argc, &argv, NULL);

  v2r_setup_test_ndsphere();
  v2r_setup_test_spheroid();
  

  return g_test_run();
}
