// solve eliptic PDE(2D)
// 
// for each level k, discretized equation is following.
// a_{ix,iy,k}*u_{ix,iy,k} + b_{ix,iy,k}*u_{ix+1,iy,k} + c_{ix,iy,k}*u_{ix-1,iy}
//                  + d_{ix,iy,k}*u_{ix,iy+1} + e_{ix,iy,k}*u_{ix,iy-1,k} = right{ix,iy,k}
//
// Assumption
//  - Dimention = 2^(max_level) * coasest_dimention + 2^(max_level) - 1
//  - Square domain
//  - dirichlet condition
//  - V-cycle multi-grid
//  - Gauss-Seidel Method for smoother

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "example.hpp"

double eps = 1e-6;
const int gauss_iter = 1;

const double L = 2.0; // the size of domain is L x L
const int coarsest_dim = 2; // the number of observing point in the domain
const int maxlevel =  9; // level = 0,1,...,maxlevel (fine -> coarse)
const int finest_dim = std::pow(2, maxlevel) * coarsest_dim + std::pow(2, maxlevel) - 1;
// coefficients for each level
std::vector<std::vector<std::vector<double>>> a(maxlevel+1),b(maxlevel+1),c(maxlevel+1),
                                 d(maxlevel+1),e(maxlevel+1),right(maxlevel+1);
                                
// solutions for each level
std::vector<std::vector<std::vector<double>>> u(maxlevel+1);

// residuals for each level
std::vector<std::vector<std::vector<double>>> residual(maxlevel+1);

// numbers for modifying solution
std::vector<std::vector<std::vector<double>>> delta_u(maxlevel + 1);


void set_initial()
{
    std::cout << "[SEQUENCE]: allocate" << std::endl;
    // allocation
    int dim = coarsest_dim;
    for(int level = maxlevel; level >= 0; --level)
    {
        printf("%d x %d\n", dim, dim);
        // int dim2D_pad = (dim + 2)*(dim + 2); // 4*dim + 4 is needed for boundary
        // coeficients
        a[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        b[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        c[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        d[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        e[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        // variants for calculation
        right[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        delta_u[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        residual[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));
        u[level].resize(dim + 2, std::vector<double>(dim + 2, 0.0));

        dim = 2 * dim + 1;
    }

    std::cout << "[SEQUENCE]: set coefficients" << std::endl;
    // set coefficients
    poisson_example(maxlevel, coarsest_dim, L, a, b, c, d, e, right);
}

// FUTURE: multi-color(RED-BLACK)
void gauss_seidel(int level, int dim)
{
    // a loop of gauss-seidel
    for(int iter = 0; iter < gauss_iter; ++iter)
    {
        for(int ix = 1; ix <= dim; ++ix)
        {
            for(int iy = 1; iy <= dim; ++iy)
            {
                u[level][ix][iy] = 1.0/a[level][ix][iy]
                                    * (
                                        right[level][ix][iy]
                                        - b[level][ix][iy]*u[level][ix+1][iy]
                                        - c[level][ix][iy]*u[level][ix-1][iy]
                                        - d[level][ix][iy]*u[level][ix][iy+1]
                                        - e[level][ix][iy]*u[level][ix][iy-1]
                                    );
            }
        }
    }

    // calculate residual
    for(int ix = 1; ix <= dim; ++ix)
    {
        for(int iy = 1; iy <= dim; ++iy)
        {
            residual[level][ix][iy] = right[level][ix][iy]
                                    - a[level][ix][iy] * u[level][ix][iy]
                                    - b[level][ix][iy] * u[level][ix+1][iy]
                                    - c[level][ix][iy] * u[level][ix-1][iy]
                                    - d[level][ix][iy] * u[level][ix][iy+1]
                                    - e[level][ix][iy] * u[level][ix][iy-1];
        }
    }
}

// 2^n * h -> 2^(n+1) h
void interp_restrict(int level, int dim)
{
    int next_dim = (dim - 1)/2;

    for(int ix = 1; ix <= next_dim; ++ix)    
    {
        for(int iy = 1; iy <= next_dim; ++iy)
        {
            right[level+1][ix][iy] = 0.125 *   (   residual[level][2*ix+1][2*iy]
                                                +  residual[level][2*ix-1][2*iy]
                                                +  residual[level][2*ix][2*iy+1]
                                                +  residual[level][2*ix][2*iy-1])
                                   + 0.25  *       residual[level][2*ix][2*iy]
                                   + 0.0625*   ( residual[level][2*ix+1][2*iy+1]
                                                +residual[level][2*ix-1][2*iy+1]
                                                +residual[level][2*ix+1][2*iy-1]
                                                +residual[level][2*ix-1][2*iy-1]);
        }
    }
}

// 2^(n+1) * h -> 2^n * h
void interp_prolong(const int level, const int dim)
{
    for(int ix = 0; ix <= dim; ++ix)
    {
        for(int iy = 0; iy <= dim; ++iy)        
        {
            delta_u[level-1][2*ix][2*iy] = u[level][ix][iy];
            delta_u[level-1][2*ix][2*iy+1] = 0.5 * (u[level][ix][iy] + u[level][ix][iy+1]);
            delta_u[level-1][2*ix+1][2*iy] = 0.5 * (u[level][ix+1][iy] + u[level][ix][iy]);
            delta_u[level-1][2*ix+1][2*iy+1] = 0.25 * ( u[level][ix][iy] + u[level][ix+1][iy]
                                                      + u[level][ix][iy+1] + u[level][ix+1][iy+1]);
        }
    }


    int next_dim = 2*dim + 1;
    for(int ix = 1; ix <= next_dim; ++ix)
    {
        for(int iy = 1; iy <= next_dim; ++iy)        
        {
            u[level-1][ix][iy] += delta_u[level-1][ix][iy];
        }
    }
}

double calculate_error()
{
    double sum_error = 0.0;

    // calculate residual
    for(int ix = 1; ix <= finest_dim; ++ix)
    {
        for(int iy = 1; iy <= finest_dim; ++iy)
        {
            double error = right[0][ix][iy]
                            - a[0][ix][iy]*u[0][ix][iy]  
                            - b[0][ix][iy]*u[0][ix+1][iy]
                            - c[0][ix][iy]*u[0][ix-1][iy]
                            - d[0][ix][iy]*u[0][ix][iy+1]
                            - e[0][ix][iy]*u[0][ix][iy-1];
            sum_error += error * error;
        }
    }

    return std::sqrt(sum_error);
}

void output_progress(int iter, double error)
{
    std::printf("\riter: %d, error:%15.11f", iter, error);
}

void fill_zero(int begin_level, int end_level, std::vector<std::vector<std::vector<double>>>& ary)
{
    for(int level = begin_level; level <= end_level; ++level)        
    {
        for(auto& column : ary[level])
        {
            for(auto& factor : column)
            {
                factor = 0.0;
            }
        }
    }
}

void show_array(const std::vector<std::vector<std::vector<double>>>& ary)
{
    std::cout << "============================"<<  std::endl;
    for(const auto& level_matrix : ary)
    {
        for(const auto& columns : level_matrix)
        {
            for(const auto& factor : columns)
            {
                std::printf("%f,", factor);
            }
            std::cout << std::endl;
        }
        std::cout << std::endl;
    }
    std::cout << "============================"<<  std::endl;
}

void debug_report()
{
    show_array(delta_u);
    std::exit(1);
}

void calculate()
{
    std::cout << "[SEQUENCE]: run multigrid" << std::endl;
    for(int iter = 1;;++iter)
    {
        // restriction
        int dim = finest_dim;
        gauss_seidel(0, dim);
        interp_restrict(0, dim);
        dim = (dim - 1)/2;
        for(int level = 1; level <= maxlevel-1; ++level)
        {
            gauss_seidel(level, dim);
            interp_restrict(level, dim);
            dim = (dim - 1)/2;
        }
        gauss_seidel(maxlevel, dim); // now, dim = 2^(coarsest_power)

        // prolongation
        for(int level = maxlevel; level >= 2; --level)
        {
            interp_prolong(level, dim);
            gauss_seidel(level-1, 2*dim+1);
            dim = 2 * dim + 1;
        }
        interp_prolong(1, dim);

        double error = calculate_error();
        if(error < eps)
        {
            // finish!
            return;
        }

        output_progress(iter, error);
        
        // fill zero for u
        fill_zero(1, maxlevel, u);
    }
}


int main()
{
    set_initial();
    calculate();
}
