#ifndef MLCOALSIM_STR_H
#define MLCOALSIM_STR_H

#include <stdarg.h>
#include <stdlib.h>

typedef char str;

str *str_init(size_t size);
str *str_alloc(str **string, size_t size);

// Alloc a new string if there is no more room to add a substring with "size" characters.
// The new string will have a (size + size/2) cap.
int str_maybe_alloc(str **s, size_t size);
void str_free(str **s);

size_t str_get_cap(str **s);
size_t str_get_len(str **s);
void str_set_len(str **s, size_t size);
void str_set_pos(str **s, size_t size, int c);

// Push a an array of characters into an existing string, indicating the length of characters to be pushed.
// If indicated length is not enough as to expand the current string, then added content will be cutted.
// This function is not intended to be used unless you know exactly what you are doing. Use `str_push_s` instead.
// A potential use case would be to make a string large enough the first time you add content to it, so next time
// you can use `str_push_s` and it wouldn't translate into a `realloc` (e.g., `str_push_l(&s, "something", 100)`)
int str_push_l(str **s, const char *d, size_t l);
// Add new content to an existing string.
int str_push_s(str **s, const char *d);

int str_vprintf(str **s, const char *f, va_list ap);
int str_printf(str **s, const char *f, ...);

#endif // MLCOALSIM_STR_H
