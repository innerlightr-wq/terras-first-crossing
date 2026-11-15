// Exact cyclic span of 163 normalized residues for seven five-one history bytes.
// Split the final sum into disjoint ranges so that all 1.37 billion residues
// do not need to reside in memory at the same time.
#include <algorithm>
#include <bit>
#include <cassert>
#include <cstdint>
#include <iostream>
#include <limits>
#include <vector>

using U = std::uint64_t;
using Wide = __uint128_t;
constexpr U modulus = 50031545098999707ULL; // 3^35
constexpr unsigned multiplicity = 163;
constexpr unsigned chunk_count = 1024;

U mul(U a, U b) { return Wide(a) * b % modulus; }

U power(U base, unsigned exponent) {
    U value = 1;
    while (exponent > 0) {
        if (exponent % 2) value = mul(value, base);
        base = mul(base, base);
        exponent /= 2;
    }
    return value;
}

std::vector<U> byte_carries() {
    std::vector<U> result;
    for (unsigned mask = 0; mask < 128; ++mask) {
        if (std::popcount(mask) != 5) continue;
        U carry = 0;
        for (unsigned j = 0; j < 7; ++j) {
            if ((mask >> j) & 1U) carry = 3 * carry + (U(1) << j);
        }
        result.push_back(carry);
    }
    assert(result.size() == 21);
    return result;
}

void unique_sort(std::vector<U>& values) {
    std::sort(values.begin(), values.end());
    values.erase(std::unique(values.begin(), values.end()), values.end());
}

struct Chunk {
    U count, span;
    std::vector<U> head, tail;
};

int main() {
    constexpr unsigned history = 7, forward_bits = 104;
    const U normalization = power((modulus + 1) / 2, 8 * history + forward_bits);
    const auto carries = byte_carries();
    std::vector<U> residues = {0}, offsets;
    for (unsigned j = 0; j < history; ++j) {
        U scale = mul(normalization, mul(power(243, history - 1 - j), power(256, j)));
        offsets.clear();
        for (U carry : carries) offsets.push_back(mul(carry, scale));
        if (j + 1 == history) break;
        std::vector<U> next;
        next.reserve(residues.size() * offsets.size());
        for (U offset : offsets)
            for (U x : residues) next.push_back((x + offset) % modulus);
        unique_sort(next);
        residues.swap(next);
    }

    std::vector<Chunk> chunks(chunk_count);
    #pragma omp parallel for schedule(dynamic)
    for (unsigned i = 0; i < chunk_count; ++i) {
        const U lo = Wide(modulus) * i / chunk_count;
        const U hi = Wide(modulus) * (i + 1) / chunk_count;
        std::vector<U> values;
        for (U offset : offsets) {
            U from = (lo + modulus - offset) % modulus;
            U to = (hi + modulus - offset) % modulus;
            auto first = std::lower_bound(residues.begin(), residues.end(), from);
            auto last = std::lower_bound(residues.begin(), residues.end(), to);
            auto append = [&](auto begin, auto end) {
                for (; begin != end; ++begin)
                    values.push_back((*begin + offset) % modulus);
            };
            if (from < to) append(first, last);
            else {
                append(first, residues.end());
                append(residues.begin(), last);
            }
        }
        unique_sort(values);
        assert(values.size() >= multiplicity);
        assert(lo <= values.front() && values.back() < hi);
        U span = modulus;
        for (U j = multiplicity - 1; j < values.size(); ++j)
            span = std::min(span, values[j] - values[j - (multiplicity - 1)]);
        chunks[i] = {values.size(), span,
                     {values.begin(), values.begin() + multiplicity - 1},
                     {values.end() - (multiplicity - 1), values.end()}};
    }

    U count = 0, span = modulus;
    for (unsigned i = 0; i < chunk_count; ++i) {
        const auto& current = chunks[i];
        const auto& next = chunks[(i + 1) % chunk_count];
        count += current.count;
        span = std::min(span, current.span);
        for (unsigned j = 0; j < multiplicity - 1; ++j) {
            U end = next.head[j] + (i + 1 == chunk_count ? modulus : 0);
            span = std::min(span, end - current.tail[j]);
        }
    }
    assert(count == 1372257936 && span == 3611678508ULL);
    std::cout << "Distinct residues: " << count << "\nMinimum cyclic span: " << span << "\n";
}
