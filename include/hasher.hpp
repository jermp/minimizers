#pragma once

#include "../external/smhasher/src/MurmurHash2.h"
#include "../external/smhasher/src/MurmurHash2.cpp"
#include "../external/smhasher/src/City.h"
#include "../external/smhasher/src/City.cpp"

namespace minimizers {

struct murmurhash2_64 {
    typedef uint64_t hash_type;
    static uint64_t hash(char const* bytes, const uint64_t num_bytes, const uint64_t seed) {
        return MurmurHash64A(bytes, num_bytes, seed);
    }
};

struct cityhash_128 {
    typedef __uint128_t hash_type;
    static __uint128_t hash(char const* bytes, const uint64_t num_bytes, const uint64_t seed) {
        auto ret = CityHash128WithSeed(bytes, num_bytes, {seed, seed});
        __uint128_t out = 0;
        out += __uint128_t(ret.first);
        out += __uint128_t(ret.second) << 64;
        return out;
    }
};

typedef murmurhash2_64 hasher64_type;
typedef cityhash_128 hasher128_type;

}  // namespace minimizers
