#include "str.h"
#include <stdarg.h>
#include <stdio.h>
#include <string.h>

#define HDR_LEN (sizeof (size_t) * 2)
#define CAP_OFF 2
#define LEN_OFF 1

inline str *str_init(size_t size) {
    return str_alloc(NULL, size);
}

str *str_alloc(str **s, size_t size) {
    // If incoming string already exists, then we `realloc` by extending it to its current "CAP" + the new size.
    size_t *n = realloc(s ? *s ? (size_t *) *s - CAP_OFF : 0 : 0, HDR_LEN + str_get_cap(s) + size + 1);

    if (!n) {
        return NULL;
    }

    n[0] = (s ? *s ? n[0] : 0 : 0) + size + 1; // Update "CAP" offset. Now we have what we had before + the new size
    n[1] = (s ? *s ? n[1] : 0 : 0);

    // It is super important to notice that we return a pointer with offset `CAP_OFF`. First two elements
    // points to the "cap" and "length" offsets. Therefore, we can quickly access to such elements by using
    // negative indices (i.e., s[-1] and s[-2]).
    return (char *) &n[CAP_OFF];
}

int str_maybe_alloc(str **s, size_t size) {
    // Check the string will overflow when new characters added to it.
    // For example, if current lenth is 8, current content is "Hello", and 5 characters will be added,
    // then we need to "extend" the current string.
    if (str_get_len(s) + size < str_get_cap(s)) {
        return 0;
    }

    // We extend the string with the number of characters we pretent to add + the half of it.
    // For example, if we want to add 10 characters, then we extend the string by 15 characters. 
    str *n = str_alloc(s, (size + (size >> 1)));

    if (!n) {
        return -1; 
    }

    // Shouldn't we free the initial "s"?: NO, because `str_alloc` is using `realloc`. Hence, `n` and `s` are the same.
    // Since `realloc` may move the allocated memory to a new location, we need to update our variable `s`.
    *s = n;
    
    return 0;
}

inline size_t str_get_cap(str **s) {
    return s ? *s ? ((size_t *) *s)[-CAP_OFF] : 0 : 0;
}

inline size_t str_get_len(str **s) {
    return s ? *s ? ((size_t *) *s)[-LEN_OFF] : 0 : 0;
}

inline void str_set_len(str **s, size_t size) {
    ((size_t *) *s)[-LEN_OFF] = size;
    str_set_pos(s, size, 0);
}

inline void str_set_pos(str **s, size_t size, int c) {
    *(*s + size) = c;
}

void str_free(str **s) {
    if (*s) {
        free((size_t *) *s - CAP_OFF);
        *s = NULL;
    }
}

int str_push_l(str **s, const char *d, size_t l) {
    // Make room in the existing string as per the indicated new length.
    if (str_maybe_alloc(s, l) < 0) {
        return -1;
    }

    // Note: new content will be cutted if indicated new lenth is not large enough.
    memcpy(*s + str_get_len(s), d, l + 1);
    str_set_len(s, str_get_len(s) + l);

    return 0;
}

int str_push_s(str **s, const char *d) {
    return str_push_l(s, d, strlen(d));
}

int str_vprintf(str **s, const char *f, va_list ap) {
    va_list ap2;
    va_copy(ap2, ap);
    int l1 = vsnprintf(NULL, 0, f, ap2);
    va_end(ap2);

    if (l1 <= 0) {
        return -1;
    }

    if (str_maybe_alloc(s, (size_t) l1) < 0) {
        return -1;
    }

    if (vsnprintf(*s + str_get_len(s), (size_t) l1 + 1, f, ap) != l1) {
        return -1;
    }

    str_set_len(s, str_get_len(s) + (size_t) l1);
    return 0;
}

int str_printf(str **s, const char *f, ...) {
    va_list ap;
    va_start(ap, f);
    int ret = str_vprintf(s, f, ap);
    va_end(ap);
    return ret;
}
