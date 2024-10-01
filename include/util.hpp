#pragma once

namespace minimizers {

namespace constants {
constexpr uint64_t default_hash_size = 128;
constexpr uint64_t default_seed = 1234567890;
}  // namespace constants

namespace util {

static void print_cmd(int argc, char** argv) {
    for (int i = 0; i != argc; ++i) std::cout << argv[i] << ' ';
    std::cout << std::endl;
}

static bool begins_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.begin());
}

static bool ends_with(std::string const& str, std::string const& pattern) {
    if (pattern.size() > str.size()) return false;
    return std::equal(pattern.begin(), pattern.end(), str.end() - pattern.size());
}

}  // namespace util

}  // namespace minimizers
