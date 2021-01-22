#include "test.h"
#include <assert.h>
#include <str.h>

static const char long_str[] = "Lorem ipsum dolor sit amet, consectetur adipiscing elit. Fusce fermentum commodo mauris non sodales. Vivamus ullamcorper consectetur urna vitae fringilla. Nulla eget libero est.";

static MunitResult
test_str_init(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_init(10);

    munit_assert_size(str_get_cap(&s), ==, 11);
    munit_assert_size(str_get_len(&s), ==, 0);
    munit_assert_not_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_alloc(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *base = str_init(10);
    str *s = str_alloc(&base, 20);

    munit_assert_size(str_get_cap(&s), ==, 32); // 32 = 10 from base + 20 from s + 2 from two extra length in each string
    munit_assert_size(str_get_len(&s), ==, 0);
    munit_assert_not_null(s);

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_maybe_alloc(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_init(10);

    munit_assert_size(str_get_cap(&s), ==, 11);

    str_maybe_alloc(&s, 10);
    munit_assert_size(str_get_cap(&s), ==, 11);

    str_maybe_alloc(&s, 20); // the string should be extended by (11 + 20 + 20/2 + 1). Extension adds 1 extra element.
    munit_assert_size(str_get_cap(&s), ==, 42);

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_alloc_null(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_alloc(NULL, 20);

    assert(str_get_cap(&s) == 21);
    assert(str_get_len(&s) == 0);
    munit_assert_not_null(s);

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_free_null(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = NULL;

    munit_assert_null(s);
    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_push_l(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_alloc(NULL, 10);

    str_push_l(&s, "hello", 5);
    munit_assert_size(str_get_cap(&s), ==, 11);
    munit_assert_size(str_get_len(&s), ==, 5);
    munit_assert_string_equal(s, "hello");

    str_push_l(&s, ", world", 7);
    munit_assert_size(str_get_cap(&s), ==, 11 + 7 + 7/2 + 1);
    munit_assert_size(str_get_len(&s), ==, 12);
    munit_assert_string_equal(s, "hello, world");

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_push_l_boundaries(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_alloc(NULL, 1);

    str_push_l(&s, "hello", 1);
    munit_assert_size(str_get_cap(&s), ==, 2);
    munit_assert_size(str_get_len(&s), ==, 1);
    munit_assert_string_equal(s, "h");

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_push_s(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = NULL;

    str_push_s(&s, long_str);
    size_t length = strlen(long_str);
    munit_assert_size(str_get_cap(&s), ==, length + length/2 + 1);
    munit_assert_size(str_get_len(&s), ==, length);
    munit_assert_string_equal(s, long_str);

    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

static MunitResult
test_str_printf(MUNIT_UNUSED const MunitParameter params[], MUNIT_UNUSED void *fixture) {
    str *s = str_init(100);

    str_push_s(&s, "Hello World!!");

    str_printf(&s, " %d", 12345);

    munit_assert_string_equal(s, "Hello World!! 12345");
    
    str_free(&s);
    munit_assert_null(s);

    return MUNIT_OK;
}

MunitTest tests_str[] = {
    test_basic(test_str_init),
    test_basic(test_str_alloc),
    test_basic(test_str_alloc_null),
    test_basic(test_str_maybe_alloc),
    test_basic(test_str_free_null),
    test_basic(test_str_push_l),
    test_basic(test_str_push_l_boundaries),
    test_basic(test_str_push_s),
    test_basic(test_str_printf),
    test_endtests,
};

MunitSuite suite_str = {
    "/strings",
    tests_str,
    NULL,
    1,
    MUNIT_SUITE_OPTION_NONE,
};
