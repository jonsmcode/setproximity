#include <unordered_set>
#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>
#include "hashfunctions.hpp"
#include "util.hpp"


static size_t set_intersection_size(const std::unordered_set<uint64_t> &set_a,
                                    const std::unordered_set<uint64_t> &set_b)
{
    if(set_b.size() < set_a.size()) {
        return set_intersection_size(set_b, set_a);
    }
    size_t cardinality = 0;
    for(std::unordered_set<uint64_t>::const_iterator it = set_a.begin(); it != set_a.end(); it++) {
        if(set_b.find(*it) != set_b.end())
            cardinality++;
    }
    return cardinality;
}

static size_t set_union_size(const std::unordered_set<uint64_t> &set_a,
                             const std::unordered_set<uint64_t> &set_b)
{
    std::unordered_set<uint64_t> set_union;
    set_union.insert(set_a.begin(), set_a.end());
    set_union.insert(set_b.begin(), set_b.end());
    return set_union.size();
}


std::vector<std::unordered_set<uint64_t>> fill_hts(const std::filesystem::path &filepath, const double percentage)
{
    std::vector<std::unordered_set<uint64_t>> kmersets;
    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});
    for(auto & record : stream) {
        std::unordered_set<uint64_t> kmerset;
        for(auto && kmer : record.sequence() | kmer_view) {
            if((double) std::rand()/RAND_MAX <= percentage)
                kmerset.insert(kmer);
        }
        kmersets.push_back(kmerset);
    }
    return kmersets;
}


double jaccard_similarity(const std::unordered_set<uint64_t> &kmerset_a,
                          const std::unordered_set<uint64_t> &kmerset_b)
{
    return (double) set_intersection_size(kmerset_a, kmerset_b)/set_union_size(kmerset_a, kmerset_b);
}


void jaccard_similarities(std::vector<std::unordered_set<uint64_t>> &kmersets, double* matrix, const size_t n)
{
    for(int i = 0; i < n; i++) {
        matrix[i*(n+1)] = 1;
        for(int j = i+1; j < n; j++) {
            matrix[i*n+j] = matrix[j*n+i] = jaccard_similarity(kmersets[i], kmersets[j]);
        }
    }
}


int main(int argc, char** argv)
{
    if(argc != 3) {
        std::cout << "usage: provide file and sample probability in [0,1]\n";
        return -1;
    }

    std::filesystem::path filepath = argv[1];
    double percentage = std::stod(argv[2]);

    std::vector<std::unordered_set<uint64_t>> kmersets = fill_hts(filepath, percentage);

    const size_t n = kmersets.size();
    double* matrix = new double[n*n];
    jaccard_similarities(kmersets, matrix, n);

    print_matrix(matrix, n);

    delete[] matrix;
}
