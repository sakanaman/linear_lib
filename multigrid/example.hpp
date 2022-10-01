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
    return exact_poisson(x, y);
}

double rhs_poisson(double x, double y)
{
    return -exact_poisson(x, y);
}


void poisson_example(int maxlevel, int coasest_power, int finest_power, 
                     double L,
                     std::vector<std::vector<double>>& a,
                     std::vector<std::vector<double>>& b,
                     std::vector<std::vector<double>>& c,
                     std::vector<std::vector<double>>& d,
                     std::vector<std::vector<double>>& e,
                     std::vector<std::vector<double>>& right)
{
    int dim = std::pow(2, finest_power);
    for(int level = 0; level <= maxlevel; ++level) 
    {
        for(int ix = 1; ix <= dim; ++ix)
        {
            for(int iy = 1; iy <= dim; ++iy)
            {
                int index_coef = (ix - 1) * dim + (iy - 1);
                double h = L/(dim + 1);
                double x = h * ix;
                double y = h * iy;

                a[level][index_coef] = -4.0/(h*h); // center
                b[level][index_coef] = 1.0/(h*h); // right
                c[level][index_coef] = 1.0/(h*h); // left
                d[level][index_coef] = 1.0/(h*h); // up
                e[level][index_coef] = 1.0/(h*h); // down
                right[level][index_coef] = rhs_poisson(x,y);
                
                if(ix == 1)
                {
                    right[level][index_coef] -= c[level][index_coef] * boundary_poisson(0,y);
                    c[level][index_coef] = 0.0;
                }
                if(iy == 1)
                {
                    right[level][index_coef] -= e[level][index_coef] * boundary_poisson(x,0);
                    e[level][index_coef] = 0.0;
                }
                if(ix == dim)
                {
                    right[level][index_coef] -= b[level][index_coef] * boundary_poisson(L,y);
                    b[level][index_coef] = 0.0;
                }
                if(iy == dim)
                {
                    right[level][index_coef] -= d[level][index_coef] * boundary_poisson(x,L);
                }
            }
        }
        dim /= 2;
    }
}

#endif