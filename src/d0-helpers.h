/**
 * \file  d0-helpers.h
 * \brief Functions that are common to all donor calculations
 */
#ifndef D0_HELPERS
#define D0_HELPERS

#include <stdlib.h>
#include "struct.h"

data11 * read_v      (size_t       *n);
double   read_delta_z(data11       *Vp);
double   V_min       (data11       *Vp,
                      const size_t  n);
#endif
