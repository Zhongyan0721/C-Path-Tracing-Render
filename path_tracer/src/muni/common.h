#pragma once
#include <spdlog/spdlog.h>
#include <linalg.h>

namespace muni {

// linalg::vec<T,M> defines a fixed-length vector containing exactly M elements of type T
template<int N, class T> using Vec = linalg::vec<T, N>;
template<class T> using Vec3 = Vec<3, T>;
template<class T> using Vec2 = Vec<2, T>;


using Vec3f = Vec3<float>;
using Vec2f = Vec2<float>;


// linalg::mat<T,M,N> defines a fixed-size matrix containing exactly M rows and N columns of type T, in column-major order.
template<class T, int M, int N> using Mat = linalg::mat<T, M, N>;
template<class T> using Mat2 = Mat<T, 2, 2>;
template<class T> using Mat3 = Mat<T, 3, 3>;
template<class T> using Mat4 = Mat<T, 4, 4>;

using Mat2f = Mat2<float>;
using Mat3f = Mat3<float>;
using Mat4f = Mat4<float>;
}  // namespace muni

// Base class for both vec and mat fmtlib formatters.
// Based on the great blog tutorial: https://wgml.pl/blog/formatting-user-defined-types-fmt.html
template<typename V, typename T, bool Newline> struct vecmat_formatter {
    using underlying_formatter_type = fmt::formatter<T>;

    template<typename ParseContext> constexpr auto parse(ParseContext &ctx) {
        return underlying_formatter.parse(ctx);
    }

    template<typename FormatContext>
    auto format(const V &v, FormatContext &ctx) {
        fmt::format_to(ctx.out(), "{{");
        auto it = begin(v);
        while (true) {
            ctx.advance_to(underlying_formatter.format(*it, ctx));
            if (++it == end(v)) {
                fmt::format_to(ctx.out(), "}}");
                break;
            } else
                fmt::format_to(ctx.out(), ",{} ", Newline ? "\n" : "");
        }
        return ctx.out();
    }

protected:
    underlying_formatter_type underlying_formatter;
};

template<typename T, int N>
struct fmt::formatter<linalg::vec<T, N>>
    : public vecmat_formatter<linalg::vec<T, N>, T, false> {};

template<typename T, int M, int N>
struct fmt::formatter<linalg::mat<T, M, N>>
    : public vecmat_formatter<linalg::mat<T, M, N>, linalg::vec<T, N>, true> {};
