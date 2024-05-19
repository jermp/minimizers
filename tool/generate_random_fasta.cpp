#include <iostream>
#include <fstream>
#include <ctime>
#include <cstdint>

#include "external/cmd_line_parser/include/parser.hpp"

int main(int argc, char** argv) {
    cmd_line_parser::parser parser(argc, argv);
    parser.add("output_filename", "Name of the output file (should end with .fa or .fasta).", "-o",
               true);
    parser.add("length", "Length of the random sequence.", "-n", true);

    if (!parser.parse()) return 1;

    const auto output_filename = parser.get<std::string>("output_filename");
    const auto n = parser.get<uint64_t>("length");

    std::ofstream out(output_filename.c_str());
    if (!out.good()) {
        std::cerr << "Error in opening file" << std::endl;
        return 1;
    }
    out << "<\n";  // fasta header

    srand(time(NULL));
    const std::string alphabet("ACGT");

    for (uint64_t i = 0; i != n; ++i) { out << alphabet[rand() % 4]; }
    out << '\n';

    out.close();

    return 0;
}
