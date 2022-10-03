#ifndef EXAMPLE_HPP
#define EXAMPLE_HPP

#include <cmath>
#include <vector>

double exact_poisson(double x, double y)
{
    return std::cos(x) + std::cos(y);
}

double boundary_poisson(double x, double y)
{
    return 1.0;
}

double rhs_poisson(double x, double y)
{
    return 1.0;
}


void poisson_example(int maxlevel, int coasest_dim, 
                     double L,
                     std::vector<std::vector<std::vector<double>>>& a,
                     std::vector<std::vector<std::vector<double>>>& b,
                     std::vector<std::vector<std::vector<double>>>& c,
                     std::vector<std::vector<std::vector<double>>>& d,
                     std::vector<std::vector<std::vector<double>>>& e,
                     std::vector<std::vector<std::vector<double>>>& right)
{
    int dim = coasest_dim;
    for(int level = maxlevel; level >= 0; --level) 
    {
        for(int ix = 1; ix <= dim; ++ix)
        {
            for(int iy = 1; iy <= dim; ++iy)
            {
                double h = L/(dim + 1);
                double x = h * ix;
                double y = h * iy;

                a[level][ix][iy] = -4.0/(h*h); // center
                b[level][ix][iy] = 1.0/(h*h); // right
                c[level][ix][iy] = 1.0/(h*h); // left
                d[level][ix][iy] = 1.0/(h*h); // up
                e[level][ix][iy] = 1.0/(h*h); // down

                if(level == 0)
                {
                    right[level][ix][iy] = rhs_poisson(x,y); // only for level = 0, but not consider.
                }

                if(ix == 1)
                {
                    if(level == 0)
                    {
                        right[level][ix][iy] -=  1.0/(h*h)* boundary_poisson(0,y); // only for level = 0, but not consider.
                    }
                    c[level][ix][iy] = 0.0;
                }
                if(iy == 1)
                {
                    if(level == 0)
                    {
                        right[level][ix][iy] -= 1.0/(h*h) * boundary_poisson(x,0);// only for level = 0, but not consider.
                    }
                    e[level][ix][iy] = 0.0;
                }
                if(ix == dim)
                {
                    if(level == 0)
                    {
                        right[level][ix][iy] -= 1.0/(h*h) * boundary_poisson(L,y);// only for level = 0, but not consider.
                    }
                    b[level][ix][iy] = 0.0;
                }
                if(iy == dim)
                {
                    if(level == 0)
                    {
                        right[level][ix][iy] -= 1.0/(h*h) * boundary_poisson(x,L);// only for level = 0, but not consider.
                    }
                    d[level][ix][iy] = 0.0;
                }
            }
        }
        dim = 2 * dim + 1;
    }
}

#endif