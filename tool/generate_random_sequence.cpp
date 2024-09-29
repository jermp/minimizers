#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdint>

#include "external/cmd_line_parser/include/parser.hpp"

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("output_filename", "Name of the output file.", "-o", true);
    parser.add("length", "Length of the random sequence.", "-n", true);
    parser.add("alphabet_size", "Size of the alphabet (must be <= 256).", "-s", true);

    if (!parser.parse()) return 1;

    const auto output_filename = parser.get<std::string>("output_filename");
    const uint64_t n = parser.get<uint64_t>("length");
    const uint64_t alphabet_size = parser.get<uint64_t>("alphabet_size");
    if (alphabet_size > 256) {
        std::cerr << "Size of the alphabet must be <= 256." << std::endl;
        return 1;
    }

    std::ofstream out(output_filename.c_str(), std::ofstream::binary);
    if (!out.good()) {
        std::cerr << "Error in opening file" << std::endl;
        return 1;
    }

    /* header: sequence length and alphabet size */
    out.write(reinterpret_cast<char const*>(&n), sizeof(n));
    out.write(reinterpret_cast<char const*>(&alphabet_size), sizeof(alphabet_size));

    srand(time(NULL));
    for (uint64_t i = 0; i != n; ++i) out << char(rand() % alphabet_size);
    out.close();

    return 0;
}
