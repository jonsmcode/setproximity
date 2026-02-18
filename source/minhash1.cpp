#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"
#include "util.hpp"


void minhash_sketches(const std::filesystem::path &filepath, std::vector<std::vector<uint64_t>> &sketches,
                      const int permutations=100, uint64_t (*hashFunc)(uint64_t)=wyhash)
{
    // std::vector<uint64_t> a(permutations);
    // std::vector<uint64_t> b(permutations);
    // for(int i = 0; i < permutations; i++) {
    //     a[i] = (1+std::rand()) % prime;
    //     b[i] = std::rand() % prime;
    // }
    murmurhash2_64 hasher;
    std::vector<uint64_t> seeds(permutations);
    for(int i = 0; i < permutations; i++)
        seeds[i] = std::rand();

    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        std::vector<uint64_t> minhash_sketch(permutations, UINT64_MAX);
        for(auto && kmer : record.sequence() | kmer_view) {
            for(int i = 0; i < permutations; i++) {
                // minhash_sketch[i] = std::min((kmer*a[i]+b[i]) % prime, minhash_sketch[i]);
                const uint64_t hash = hasher.hash(kmer, seeds[i]);
                minhash_sketch[i] = std::min(hash, minhash_sketch[i]);
            }
        }
        sketches.push_back(minhash_sketch);
    }
}


double minhash_similarity(const std::vector<uint64_t> &minhashs_a,
                          const std::vector<uint64_t> &minhashs_b)
{
    const size_t permutations = minhashs_a.size();
    int y = 0;
    for(size_t i = 0; i < permutations; i++)
        y += minhashs_a[i] == minhashs_b[i];

    return (double) y/permutations;
}


void minhash_similarities(const std::vector<std::vector<uint64_t>> &sketches, double* matrix, const int n)
{
    for(int i = 0; i < n; i++) {
        matrix[i*(n+1)] = 1.0;
        for(int j = i+1; j < n; j++) {
            matrix[i*n+j] = matrix[j*n+i] = minhash_similarity(sketches[i], sketches[j]);
        }
    }

}


int main(int argc, char** argv)
{
    std::vector<std::vector<uint64_t>> sketches;
    std::filesystem::path filepath = argv[1];

    if(argc == 2) {
        minhash_sketches(filepath, sketches);
    }
    else if(argc == 3) {
        const int permutations = std::stoi(argv[2]);
        minhash_sketches(filepath, sketches, permutations);
    }
    else if(argc == 4) {
        const int permutations = std::stoi(argv[1]);
        if(std::string(argv[3]) == "stdhash")
            minhash_sketches(filepath, sketches, permutations, stdhash);
        else if(std::string(argv[3]) == "wyhash")
            minhash_sketches(filepath, sketches, permutations, wyhash);
        else if(std::string(argv[3]) == "farmhash")
            minhash_sketches(filepath, sketches, permutations, farmhash);
        else
            std::cout << "no such hash function\n";
    }
    else {
        std::cout << "usage: optionally provide number of permutations and hash function\n";
        return -1;
    }

    const size_t n = sketches.size();
    double* matrix = new double[n*n];

    minhash_similarities(sketches, matrix, n);
    
    print_matrix(matrix, n);

}
