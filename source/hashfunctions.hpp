#include "farmhash.h"

static inline uint64_t stdhash(const uint64_t x) {
    std::hash<uint64_t> hasher;
    return hasher(x);
}

static inline uint64_t wyhash(const uint64_t x) {
    __uint128_t result = x;
    result *= 0x9E3779B97F4A7C15ULL;
    return static_cast<uint64_t>(result) ^ static_cast<uint64_t>(result >> 64);
}

static inline uint64_t farmhash(uint64_t x) {
    return util::Fingerprint(x);
}



static uint64_t MurmurHash2_64(void const* key, size_t len, uint64_t seed) {
    const uint64_t m = 0xc6a4a7935bd1e995ULL;
    const int r = 47;

    uint64_t h = seed ^ (len * m);

#if defined(__arm) || defined(__arm__)
    const size_t ksize = sizeof(uint64_t);
    const unsigned char* data = (const unsigned char*)key;
    const unsigned char* end = data + (std::size_t)(len / 8) * ksize;
#else
    const uint64_t* data = (const uint64_t*)key;
    const uint64_t* end = data + (len / 8);
#endif

    while (data != end) {
#if defined(__arm) || defined(__arm__)
        uint64_t k;
        memcpy(&k, data, ksize);
        data += ksize;
#else
        uint64_t k = *data++;
#endif

        k *= m;
        k ^= k >> r;
        k *= m;

        h ^= k;
        h *= m;
    }

    const unsigned char* data2 = (const unsigned char*)data;

    switch (len & 7) {
        // fall through
        case 7:
            h ^= uint64_t(data2[6]) << 48;
        // fall through
        case 6:
            h ^= uint64_t(data2[5]) << 40;
        // fall through
        case 5:
            h ^= uint64_t(data2[4]) << 32;
        // fall through
        case 4:
            h ^= uint64_t(data2[3]) << 24;
        // fall through
        case 3:
            h ^= uint64_t(data2[2]) << 16;
        // fall through
        case 2:
            h ^= uint64_t(data2[1]) << 8;
        // fall through
        case 1:
            h ^= uint64_t(data2[0]);
            h *= m;
    };

    h ^= h >> r;
    h *= m;
    h ^= h >> r;

    return h;
}

struct murmurhash2_64 {
    // generic range of bytes
    // static inline uint64_t hash(byte_range range, uint64_t seed) {
    //     return MurmurHash2_64(range.begin, range.end - range.begin, seed);
    // }

    // specialization for std::string
    static inline uint64_t hash(std::string const& val, uint64_t seed) {
        return MurmurHash2_64(val.data(), val.size(), seed);
    }

    // specialization for uint64_t
    static inline uint64_t hash(uint64_t val, uint64_t seed) {
        return MurmurHash2_64(reinterpret_cast<char const*>(&val), sizeof(val), seed);
    }
};