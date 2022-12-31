#include "surface.h"
#include "photonmap.h"
#include "render.h"
#include "kdtree.h"

#include <vector>
#include <complex>
#include <functional>
#include <algorithm>
#include <iostream>

/* 
 * 1. Load object surfaces as points
 * 2. Create trees for objects and light sources
 * 3. Find nearest points on lightsource tree from object tree
 * 4. Project light source onto objects to define shadows
 * 5. Increase point density around shadows
 * 6. Calculate luminance at points
 * 7. Interpolate luminance between points
 * 8. Find nearest points on camera tree from object tree
 * 9. Back calculate luminance on camera tree from points
 * 10. Map wavelength power to RGB
 */

int main(void)
{
    Surface p;
    std::vector<data<Surface>> points({
        {p, {0, 1, 0}},{p, {1, 0, 0}},{p, {0, 0, 0}},
        {p, {1, 2, 4}},{p, {0, 1, 6}},{p, {0, 1, 2}}
    });
    KDTree<Surface> tree(points);

    data<Surface> d = {p, {0, 1, 3}};

    tree.insert(d);

    data<Surface> n = tree.nearestNeighbor({0, 1, 5});

    for (int i = 0; i < (int)n.p.size(); i++) {
        std::cout << n.p[i] << ' ';
    }

    std::cout << std::endl << sizeof(Surface);

    return 0;
}