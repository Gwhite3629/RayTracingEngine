#ifndef _RENDER_H_
#define _RENDER_H_

#include <vector>

// Render equation, calculates light intensity from a source at a destination
template<class T>
T Intensity(std::vector<T> destination, std::vector<T> source)
{
    return g(destination, source) * (emit(destination, source) + scattered(destination, source));
}

#endif // _RENDER_H_