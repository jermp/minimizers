#include "algorithms.hpp"

namespace minimizers {

namespace detail {

static globals_t globals;
globals_t& get_globals() { return globals; }
}  // namespace detail

double redundancy_in_density_in_perc(const double density, const double lower_bound) {
    return (density / lower_bound - 1) * 100.0;
}

double redundancy_in_density_as_factor(const double density, const double lower_bound) {
    return density / lower_bound;
}

bool is_not_forward(const uint64_t k, const uint64_t w, const uint64_t t) {
    assert(w >= 2);
    assert(t <= k);
    /*
        We know a scheme is *not* foward when there exist x and y
        such that x mod w + 1 < y mod w, where x and y are the
        positions of the smallest t-mer in window i and i-1 respectively,
        for some i > 0. So we derive: x mod w < w - 2.
        Since x is in [0..l - t] = [0..w + k - 1 - t], then x is at most
        w + k - t - 1, i.e., w + k - t - 1 mod w < w - 2.

        All possible backward jumps (y mod w, x mod w), of length y-x-1,
        are for y in [x+1..w-1].

        Note: in math, we would write (k - t - 1) mod w < w - 2,
        but here we always sum w to avoid having to take the
        modulo of negative integers when t = k.
    */
    return (w + k - t - 1) % w < w - 2;
}

/* This ignores (asymptotic) lower order terms. */
double closed_form_density(std::string const& scheme_name,                        //
                           const uint64_t k, const uint64_t w, const uint64_t t)  //
{
    if (scheme_name == "miniception") {
        return 1.67 / w;
    } else if (scheme_name == "mod_sampling") {
        bool ok = (w + k - 1 - t) % w == w - 1;
        double correction = ok ? 0 : floor(1.0 + double(k - 1.0 - t) / w) / (w + k - t);
        return double(floor(1.0 + double(k - t - 1.0) / w) + 2.0 - correction) / (w + k - t + 1.0);
    } else {
        throw std::runtime_error("unknown scheme name");
    }
}

bool is_decycling_original(char const* kmer, const uint64_t k) {
    long double im = 0;
    for (uint64_t i = 0; i != k; ++i) { im += detail::globals.sines[i] * kmer[i]; }
    long double im_rot = 0;
    for (uint64_t i = 0; i != k; ++i) { im_rot += detail::globals.sines[i + 1] * kmer[i]; }
    // std::cerr << "Im: " << im << " im_rot: " << im_rot <<
    // std::endl;
    return im > 0 and im_rot <= 0;
}

bool is_decycling_arg_pos(char const* kmer, const uint64_t k) {
    std::complex<long double> x = 0;
    for (uint64_t i = 0; i != k; ++i) { x += detail::globals.roots[i] * (long double)kmer[i]; }
    auto a = std::arg(x);
    // std::cerr << "arg:    " << a << "\nthresh: " << pi - 2 * pi / k << endl;
    return detail::globals.pi - 2 * detail::globals.pi / k < a;
}

bool is_decycling_arg_neg(char const* kmer, const uint64_t k) {
    std::complex<long double> x = 0;
    for (uint64_t i = 0; i != k; ++i) { x += detail::globals.roots[i] * (long double)kmer[i]; }
    auto a = std::arg(x);
    return -2 * detail::globals.pi / k < a and a <= 0;
}

}  // namespace minimizers
