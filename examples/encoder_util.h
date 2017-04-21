// Utility functions used by encoder binaries.

#ifndef EXAMPLES_ENCODER_UTIL_H_
#define EXAMPLES_ENCODER_UTIL_H_

#include "./aom_config.h"
#include "aom/aom_image.h"

// Returns mismatch location (?loc[0],?loc[1]) and the values at that location
// in img1 (?loc[2]) and img2 (?loc[3]).
#if CONFIG_HIGHBITDEPTH
void aom_find_mismatch_high(const aom_image_t *const img1,
                            const aom_image_t *const img2, int yloc[4],
                            int uloc[4], int vloc[4]);
#endif  // CONFIG_HIGHBITDEPTH

void aom_find_mismatch(const aom_image_t *const img1,
                       const aom_image_t *const img2, int yloc[4], int uloc[4],
                       int vloc[4]);

// Returns 1 if the two images match.
int aom_compare_img(const aom_image_t *const img1,
                    const aom_image_t *const img2);

#endif  // EXAMPLES_ENCODER_UTIL_H_
