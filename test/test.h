#ifndef MLCOALSIM_TEST_H__
#define MLCOALSIM_TEST_H__

#include <munit.h>

extern MunitSuite suite_str;

// Setup function to be passed to test functions
void *test_setup(const MunitParameter params[], void *fixture);

// Convenient macro for building a basic`MunitTest struct` with the test name, the test function, and a setup function.
#define test_basic(name)                                                        \
    { ("/" #name), name, test_setup, NULL, MUNIT_TEST_OPTION_NONE, NULL, }


// Convenient macro for building a `MunitTest struct` to include a tear-down function
#define test_endtests                                                           \
    { NULL, NULL, NULL, NULL, MUNIT_TEST_OPTION_NONE, NULL }

// Convenient macro to check a given input is a `str`
// #define test_is_str(c)                                                       \
//     do {                                                                       \
//         grim_object z = (c);                                                   \
//         munit_assert_int(grim_type(z), ==, GRIM_STRING);                       \
//         munit_assert_int(grim_direct_tag(z), ==, GRIM_INDIRECT_TAG);           \
//         munit_assert_int(I_tag(z), ==, GRIM_STRING_TAG);                       \
//     } while (0)

// #define test_check_str(c, l, v)                                              \
//     do {                                                                       \
//         grim_object y = (c);                                                   \
//         test_is_str(y);                                                      \
//         size_t L = (l);                                                        \
//         uint8_t *w = (uint8_t *)(v);                                           \
//         munit_assert_ullong(I_strlen(y), ==, L);                               \
//         munit_assert_ullong(I_strlen(y), ==, L);                               \
//         for (size_t i = 0; i < L; i++)                                         \
//             munit_assert_uint8(I_str(y)[i], ==, w[i]);                         \
//     } while (0)

#endif /* MLCOALSIM_TEST_H__ */
