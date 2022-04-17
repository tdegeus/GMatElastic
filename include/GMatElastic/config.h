/**
\file
\copyright Copyright. Tom de Geus. All rights reserved.
\license This project is released under the MIT License.
*/

#ifndef GMATELASTIC_H
#define GMATELASTIC_H

/**
\cond
*/
#define Q(x) #x
#define QUOTE(x) Q(x)

#define GMATELASTIC_ASSERT_IMPL(expr, file, line) \
    if (!(expr)) { \
        throw std::runtime_error( \
            std::string(file) + ':' + std::to_string(line) + \
            ": assertion failed (" #expr ") \n\t"); \
    }


/**
\endcond
*/

/**
All assertions are implementation as::

    GMATELASTIC_ASSERT(...)

They can be enabled by::

    #define GMATELASTIC_ENABLE_ASSERT

(before including GMatElastic).
The advantage is that:

-   File and line-number are displayed if the assertion fails.
-   GMatElastic's assertions can be enabled/disabled
    independently from those of other libraries.

\throw std::runtime_error
*/
#ifdef GMATELASTIC_ENABLE_ASSERT
#define GMATELASTIC_ASSERT(expr) GMATELASTIC_ASSERT_IMPL(expr, __FILE__, __LINE__)
#else
#define GMATELASTIC_ASSERT(expr)
#endif

/**
Linear elastic material model.
*/
namespace GMatElastic {
}

#endif
