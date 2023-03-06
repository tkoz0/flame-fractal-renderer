// define macros for system endianness

#pragma once

#if __BYTE_ORDER__ == __ORDER_LITTLE_ENDIAN__
#define TKOZ_BIG_ENDIAN 0
#define TKOZ_LITTLE_ENDIAN 1
#elif __BYTE_ORDER__ == __ORDER_BIG_ENDIAN__
#define TKOZ_BIG_ENDIAN 1
#define TKOZ_LITTLE_ENDIAN 0
#else
#error "unknown endian"
#endif
