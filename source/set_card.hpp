#include <seqan3/io/sequence_file/all.hpp>
#include <seqan3/search/views/kmer_hash.hpp>

std::vector<uint64_t> exact_sizes(const std::filesystem::path &filepath)
{
	std::vector<uint64_t> sizes;

    auto stream = seqan3::sequence_file_input<my_traits>{filepath};
    auto kmer_view = seqan3::views::kmer_hash(seqan3::ungapped{k});

    for(auto & record : stream) {
    	std::unordered_set<uint64_t> kmerset;
    	for(auto && kmer : record.sequence() | kmer_view) {
            kmerset.insert(kmer);
    	}
    	sizes.push_back(kmerset.size());
    }

	return sizes;
}