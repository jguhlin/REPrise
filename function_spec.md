# REPrise Function Specifications

This document lists all function signatures from the REPrise C++ implementation along with their Rust equivalents where available.

| Function | Return Type | Parameters | Description | Location | Rust Equivalent |
|----------|-------------|------------|-------------|----------|-----------------|
| store_cache | void | int edit_distance, vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable | Builds a cache table for all possible kmers of a given length with a specified edit distance | REPrise.hpp:6, REPrise.cpp:171 | Yes (simplified) |
| build_sortedkmers | void | priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &dist0cachetable, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable | Builds sorted kmers together with their frequencies | REPrise.hpp:7, REPrise.cpp:206 | Yes (simplified) |
| build_repeat_families | void | priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable | Builds repeat families from sorted kmers | REPrise.hpp:8, REPrise.cpp:345 | Yes (simplified) |
| findkmer | vector<seq_type> | const vector<char>& query, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable | Finds occurrences of a kmer in the sequence using the cache table | REPrise.hpp:10, REPrise.cpp:587 | Yes (simplified) |
| SA_search | void | const vector<char> &query, seq_type begin, seq_type end, char query_num, char seq_num, char rem_dist, set<tuple<seq_type, seq_type, char, char>> &matched | Performs a recursive search on the suffix array allowing mismatches | REPrise.hpp:11, REPrise.cpp:612 | Yes (simplified) |
| find_bestseed | pair<int, vector<char>> | priority_queue<pair<int, vector<char>>> &kmers, const vector<vector<tuple<seq_type, seq_type, char, char>>> &cachetable, const vector<bool> &mask_flag | Finds the best seed from the priority queue | REPrise.hpp:12, REPrise.cpp:654 | No |
| extend | int | bool isright, int seedfreq, vector<seq_type> &seed_ext | Extends consensus sequence in one direction | REPrise.hpp:13, REPrise.cpp:700 | Yes (placeholder) |
| compute_score | int | bool isright, int ext, int se, char base, const vector<vector<int>>& score_m, const vector<vector<int>>& score_ins, const vector<vector<int>>& score_del, vector<vector<vector<int>>>& score_m_bybase, vector<vector<vector<int>>>& score_ins_bybase, vector<vector<vector<int>>>& score_del_bybase | Computes alignment scores for extension | REPrise.hpp:14, REPrise.cpp:775 | Yes (simplified) |
| removetandem | void | vector<seq_type> &occs | Removes tandem occurrences closer than TANDEMDIST | REPrise.hpp:16, REPrise.cpp:510 | Yes |
| removemasked | void | vector<seq_type> &occs, const vector<bool> &mask_flag, bool isrc | Removes occurrences that overlap masked regions | REPrise.hpp:17, REPrise.cpp:525 | Yes |
| maskbyseed | void | const vector<seq_type> &occs, vector<bool> &mask_flag, bool isrc | Masks a set of occurrences | REPrise.hpp:18, REPrise.cpp:555 | Yes |
| maskbyrepeat | void | int seedfreq, const vector<seq_type> &repeatstart, const vector<seq_type> &repeatend, vector<bool> &mask_flag | Masks repeat regions | REPrise.hpp:19, REPrise.cpp:569 | No |
| maskbyrepeat_element | void | int i, seq_type elementstart, seq_type elementend, vector<bool> &mask_flag | Masks a repeat element | REPrise.hpp:20, REPrise.cpp:582 | No |
| masking_align | pair<seq_type, seq_type> | int i,seq_type consensusstart, seq_type consensusend | Performs pairwise realignment between consensus and candidates | REPrise.hpp:23, REPrise.cpp:836 | No |
| mask_extention_score | int | bool isright, int ext, int i, vector<int> &mask_score, vector<int> &mask_score_m, vector<int> &mask_score_ins, vector<int> &mask_score_del | Computes extension scores for masking alignment | REPrise.hpp:24, REPrise.cpp:919 | No |
| build_sequence | void |  | Reads sequence file and builds the padded numeric sequence | REPrise.hpp:26, REPrise.cpp:988 | Yes |
| chrtracer | pair<string, seq_type> | const seq_type stringpos | Traces chromosome information for a position | REPrise.hpp:27, REPrise.cpp:1072 | No |
| allocate_space | void |  | Allocates memory for consensus sequence | REPrise.hpp:28, REPrise.cpp:1083 | No |
| freespace | void |  | Frees allocated memory | REPrise.hpp:29, REPrise.cpp:1092 | No |
| num_to_char | char | char z | Converts numeric code to character | REPrise.hpp:30, REPrise.cpp:1098 | No |
| char_to_num | char | char c | Converts character to numeric code | REPrise.hpp:31, REPrise.cpp:1110 | Yes (private) |
| complement | char | char c | Returns complement of a nucleotide | REPrise.hpp:32, REPrise.cpp:1129 | No |
| reverse_complement | vector<char> | const vector<char> &query | Returns reverse complement of a sequence | REPrise.hpp:33, REPrise.cpp:1141 | No |
| compute_entropy | double | const vector<char> &kmer | Computes entropy of a kmer | REPrise.hpp:34, REPrise.cpp:1148 | No |
| default_k | int | seq_type len,int KMERDIST | Computes default kmer length | REPrise.hpp:35, REPrise.cpp:1170 | No |
| display_time | void | string msg | Displays timing information | REPrise.hpp:37, REPrise.cpp:1196 | No |
| print_usage | void |  | Prints usage information | REPrise.hpp:38, REPrise.cpp:1214 | No |