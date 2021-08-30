/*

(c - MIT) T.W.J. de Geus (Tom) | www.geus.me | github.com/tdegeus/GMatElastic

*/

#ifndef GMATELASTIC_H
#define GMATELASTIC_H

#ifdef GMATELASTIC_ENABLE_ASSERT

    #define GMATELASTIC_ASSERT(expr) GMATELASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
    #define GMATELASTIC_ASSERT_IMPL(expr, file, line) \
        if (!(expr)) { \
            throw std::runtime_error( \
                std::string(file) + ':' + std::to_string(line) + \
                ": assertion failed (" #expr ") \n\t"); \
        }

#else

    #define GMATELASTIC_ASSERT(expr)

#endif

#endif
