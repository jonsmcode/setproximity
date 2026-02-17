#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"

const uint8_t k = 31;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


static inline constexpr uint8_t leading_zeros(const uint64_t x, const uint8_t precision) {
    const uint64_t mask = (1ULL << precision) - 1u;
    const uint8_t rank = std::countl_zero((x << precision) | mask);
    return rank;
}

static inline constexpr uint64_t register_index(const uint64_t x, const uint8_t precision) {
    return x >> (64 - precision);
}


uint64_t hyperloglog(const std::filesystem::path &filepath, const uint8_t precision=8u, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    // TODO: implement HyperLogLog here
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    uint64_t hash;
    std::vector<uint8_t> ls;
    uint8_t bucket, temp;
    for (size_t i = 0; i < 1 << precision; i++){
        ls.push_back(0);
    }
    for(auto & record : stream) {
        for(auto && kmer : record.sequence() | kmer_view) {
            hash = hashFunc(kmer);
            temp = leading_zeros(hash, precision);
            bucket = register_index(hash, precision);
            temp = std::ranges::max(temp, ls[bucket]);
            ls[bucket] = temp;
        }
    }
    double hmean = 0;
    for (size_t i = 0; i < 1 << precision; i++){
        hmean += std::pow(2, -ls[i]);
    }
    return 0.7*(2<<precision)/hmean;
}


