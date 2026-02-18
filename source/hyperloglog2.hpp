#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"


static inline constexpr uint8_t leading_zeros(const uint64_t x, const uint8_t precision) {
    const uint64_t mask = (1ULL << precision) - 1u;
    const uint8_t rank = std::countl_zero((x << precision) | mask);
    return rank;
}

static inline constexpr uint64_t register_index(const uint64_t x, const uint8_t precision) {
    return x >> (64 - precision);
}


std::vector<uint64_t> hyperloglogs(const std::filesystem::path &filepath,
    const uint8_t precision=8u, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    std::vector<uint64_t> estimates;

    const uint64_t m = (1ULL<<precision);
    // const double alpha = 0.7213 / (1.0 + 1.079 / m);
    const double alpha = 0.7;

    uint8_t* registers = new uint8_t[m];
    std::memset(registers, 0, m*sizeof(uint8_t));

    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    for(auto & record : stream) {
        std::memset(registers, 0, m*sizeof(uint8_t));

        for(auto && kmer : record.sequence() | kmer_view) {
            const uint64_t hash = hashFunc(kmer);
            const uint64_t index = register_index(hash, precision);
            const uint8_t zeros = leading_zeros(hash, precision) + 1;
            registers[index] = std::max(registers[index], zeros);
        }

        double sum = 0.0;
        for(uint64_t j=0; j < m; ++j)
            sum += 1.0 / (1ULL << registers[j]);

        estimates.push_back(alpha*m*m/sum);
    }

    delete[] registers;

    return estimates;
}


