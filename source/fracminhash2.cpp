#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "util.hpp"
#include "hyperloglog2.hpp"
#include "set_card.hpp"


static double jaccard(const double containment, const uint64_t size_a, const uint64_t size_b) {
    return containment*size_a/(size_a+size_b-containment*size_a);
}

std::vector<std::vector<uint64_t>> fracmin_sketches(const std::filesystem::path &filepath,
    const double s, uint64_t (*hashFunc)(uint64_t)=stdhash)
{
    std::vector<std::vector<uint64_t>> sketches;
    auto fin = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : fin) {
        std::vector<uint64_t> frac_sketch;
        for(auto && kmer : record.sequence() | kmer_view) {
            uint64_t hash = hashFunc(kmer);
            if(hash <= UINT64_MAX*s)
                frac_sketch.push_back(hash);
        }
        std::sort(frac_sketch.begin(), frac_sketch.end());
        auto it = std::unique(frac_sketch.begin(), frac_sketch.end());
        frac_sketch.erase(it, frac_sketch.end());
        sketches.push_back(frac_sketch);
    }

    return sketches;
}


double fracMinHash(const std::vector<uint64_t> &frac_sketch_a, const std::vector<uint64_t> &frac_sketch_b, const uint64_t size_a,
                   const double s)
{
    std::vector<uint64_t> intersection;
    std::set_intersection(frac_sketch_a.begin(), frac_sketch_a.end(), frac_sketch_b.begin(),
                          frac_sketch_b.end(), back_inserter(intersection));
    
    return (double) intersection.size()/(frac_sketch_a.size()*(1.0-std::pow(1.0-s, size_a)));
}


void fracminhash_similarities(const std::vector<std::vector<uint64_t>> &sketches,
    const std::vector<uint64_t> &sizes, double* matrix, const int n, double s)
{
    for(int i = 0; i < n; i++) {
        matrix[i*(n+1)] = 1.0;
        for(int j = i+1; j < n; j++) {
            const double containment = fracMinHash(sketches[i], sketches[j], sizes[i], s);
            matrix[i*n+j] = matrix[j*n+i] = jaccard(containment, sizes[i], sizes[j]);
        }
    }
}


int main(int argc, char** argv)
{
    std::filesystem::path filepath = argv[1];
    double s = 0.1;
    if(argc == 3)
        s = std::stod(argv[2]);

    std::vector<std::vector<uint64_t>> sketches = fracmin_sketches(filepath, s);
    // std::vector<uint64_t> sizes = hyperloglogs(filepath, 10u);
    std::vector<uint64_t> sizes = exact_sizes(filepath);

    const size_t n = sketches.size();

    double* matrix = new double[n*n];
    fracminhash_similarities(sketches, sizes, matrix, n, s);

    print_matrix(matrix, n);

    delete[] matrix;
}
