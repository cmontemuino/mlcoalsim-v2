#include "test.h"

int main(int argc, char * const *argv) {
    MunitSuite suites[] = {
        suite_str,
    };

    MunitSuite suite = {"/mlcoalsim-str", NULL, suites, MUNIT_TEST_OPTION_SINGLE_ITERATION, MUNIT_SUITE_OPTION_NONE};
    return munit_suite_main(&suite, NULL, argc, argv);
}

void *test_setup(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    return NULL;
}
