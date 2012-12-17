#ifndef MAFIA_UTILS_H_
#define MAFIA_UTILS_H_

/** @file utils.h some utility functions and definitions */

/** compiled with device support */ 
#define MAFIA_USE_DEVICE

/** division with rounding upwards */
inline int divup(int a, int b) { return a / b + (a % b ? 1 : 0); }

/** accessing points; variables ps, n and d need be defined */
//#define PS(i, idim) ps[i * d + idim]
#define PS(i, idim) ps[idim * n + i]

#endif
