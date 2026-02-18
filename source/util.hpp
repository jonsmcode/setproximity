const uint8_t k = 15;
const uint64_t prime = (1ULL << 61) - 1u;
const int seed = 1;

struct my_traits:seqan3::sequence_file_input_default_traits_dna {
    using sequence_alphabet = seqan3::dna4;
};


void print_matrix(const double* matrix, const size_t n) {
    for(size_t i = 0; i < n; i++) {
        for(size_t j = 0; j < n; j++) {
            std::cout << matrix[i*n+j] << " ";
        }
        std::cout << '\n';
    }
}