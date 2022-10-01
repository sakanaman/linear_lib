// solve eliptic PDE(2D)
// 
// for each level k, discretized equation is following.
// a_{ix,iy,k}*u_{ix,iy,k} + b_{ix,iy,k}*u_{ix+1,iy,k} + c_{ix,iy,k}*u_{ix-1,iy}
//                  + d_{ix,iy,k}*u_{ix,iy+1} + e_{ix,iy,k}*u_{ix,iy-1,k} = right{ix,iy,k}
//
// Assumption
//  - Dimention = 2^N
//  - Square domain
//  - V-cycle multi-grid
//  - Gauss-Seidel Method for smoother

#include <vector>
#include <algorithm>
#include <cmath>
#include <iostream>
#include <string>
#include "example.hpp"

double eps = 1e-10;

const double L = 10.0; // the size of domain is L x L
const int finest_power = 9;
const int coarsest_power = 4;
const int maxlevel =  finest_power - coarsest_power; // level = 0,1,...,maxlevel
// coefficients for each level
std::vector<std::vector<double>> a(maxlevel+1),b(maxlevel+1),c(maxlevel+1),
                                 d(maxlevel+1),e(maxlevel+1),right(maxlevel+1);
                                
// solutions for each level
std::vector<std::vector<double>> u(maxlevel+1);

// residuals for each level
std::vector<std::vector<double>> residual(maxlevel+1);

// numbers for modifying solution
std::vector<std::vector<double>> delta_u(maxlevel + 1);


void set_initial()
{
    std::cout << "[SEQUENCE]: allocate" << std::endl;
    // allocation
    int dim = std::pow(2, finest_power);
    for(int power = finest_power; power >= coarsest_power; --power)
    {
        printf("%d x %d\n", dim, dim);
        int level = finest_power - power;
        int dim2D_pad = (dim + 2)*(dim + 2); // 4*dim + 4 is needed for boundary
        a[level].resize(dim2D_pad, 0.0);
        b[level].resize(dim2D_pad, 0.0);
        c[level].resize(dim2D_pad, 0.0);
        d[level].resize(dim2D_pad, 0.0);
        e[level].resize(dim2D_pad, 0.0);
        right[level].resize(dim2D_pad, 0.0);
        delta_u[level].resize(dim2D_pad, 0.0);

        // Note that these variants have different index!!
        residual[level].resize(dim2D_pad, 0.0);
        u[level].resize(dim2D_pad, 0.0);

        dim /= 2;
    }

    std::cout << "[SEQUENCE]: set coefficients" << std::endl;
    // set coefficients
    poisson_example(maxlevel, coarsest_power, finest_power, L,
                    a, b, c, d, e, right);
}

// FUTURE: multi-color(RED-BLACK)
void gauss_seidel(int level, int dim)
{
    // a loop of gauss-seidel
    for(int iter = 0; iter < 10; ++iter)
    {
        for(int ix = 1; ix <= dim; ++ix)
        {
            for(int iy = 1; iy <= dim; ++iy)
            {
                int index = ix * (dim + 2) + iy;
                int index_up = index + 1;
                int index_down = index - 1;
                int index_right = index + (dim + 2);
                int index_left = index - (dim + 2);

                u[level][index] = 1.0/a[level][index] 
                                    * (
                                        right[level][index]
                                        - b[level][index]*u[level][index_right]
                                        - c[level][index]*u[level][index_left]
                                        - d[level][index]*u[level][index_up]
                                        - e[level][index]*u[level][index_down]
                                    );
            }
        }
    }

    // calculate residual
    for(int ix = 1; ix <= dim; ++ix)
    {
        for(int iy = 1; iy <= dim; ++iy)
        {
            int index = ix * (dim + 2) + iy;
            int index_up = index + 1;
            int index_down = index - 1;
            int index_right = index + (dim + 2);
            int index_left = index - (dim + 2);

            residual[level][index] = right[level][index]
                                        - a[level][index]*u[level][index]  
                                        - b[level][index]*u[level][index_right]
                                        - c[level][index]*u[level][index_left]
                                        - d[level][index]*u[level][index_up]
                                        - e[level][index]*u[level][index_down];
        }
    }
}

void interp_restrict(int level, int dim)
{
    int next_dim = dim/2;

    for(int ix = 1; ix <= next_dim; ++ix)    
    {
        for(int iy = 1; iy <= next_dim; ++iy)
        {
            int index_coarse = ix * (next_dim + 2) + iy;

            int index_fine = 2 * ix * (dim + 2) + 2 * iy;
            int index_right = index_fine + (dim + 2);
            int index_left = index_fine - (dim + 2);
            int index_up = index_fine + 1;
            int index_down = index_fine - 1;

            right[level+1][index_coarse] = 0.125 * (residual[level][index_right] 
                                                 +  residual[level][index_left]
                                                 +  residual[level][index_up]
                                                 +  residual[level][index_down])
                                           + 0.5 *  residual[level][index_fine];
        }
    }
}

void interp_prolong(const int level, const int dim)
{
    int next_dim = dim * 2;
    for(int ix = 1; ix <= dim; ++ix)
    {
        for(int iy = 1; iy <= dim; ++iy)        
        {
            int index_fine = 2*ix * (next_dim + 2) + 2*iy;
            int index_fine_right = index_fine + (next_dim + 2);
            int index_fine_up = index_fine + 1;
            int index_fine_rightup = index_fine + (next_dim + 2) + 1;

            int index_coarse = ix * (dim + 2) + iy;
            int index_coarse_right = index_coarse + (dim + 2);
            int index_coarse_up = index_coarse + 1;
            int index_coarse_rightup = index_coarse + (dim + 2) + 1;
            
            delta_u[level-1][index_fine] = u[level][index_coarse];
            delta_u[level-1][index_fine_right] = 0.5 * (u[level][index_coarse] + u[level][index_coarse_right]);
            delta_u[level-1][index_fine_up] = 0.5 * (u[level][index_coarse] + u[level][index_coarse_up]);
            delta_u[level-1][index_fine_rightup] = 0.25 * (u[level][index_coarse] + u[level][index_coarse_right] 
                                                          +u[level][index_coarse_up] + u[level][index_coarse_rightup]);
        }
    }


    for(int ix = 1; ix <= next_dim; ++ix)
    {
        for(int iy = 1; iy <= next_dim; ++iy)        
        {
            int index_fine = ix * (next_dim + 2) + iy;

            u[level-1][index_fine] += delta_u[level-1][index_fine];
        }
    }
}

double calculate_error()
{
    double sum_error = 0.0;
    int dim = std::pow(2, finest_power);
    int level = 0;

    // calculate residual
    for(int ix = 1; ix <= dim; ++ix)
    {
        for(int iy = 1; iy <= dim; ++iy)
        {
            int index = ix * (dim + 2) + iy;
            int index_up = index + 1;
            int index_down = index - 1;
            int index_right = index + (dim + 2);
            int index_left = index - (dim + 2);

            double error = right[level][index]
                            - a[level][index]*u[level][index]  
                            - b[level][index]*u[level][index_right]
                            - c[level][index]*u[level][index_left]
                            - d[level][index]*u[level][index_up]
                            - e[level][index]*u[level][index_down];
            sum_error += error * error;
        }
    }

    return std::sqrt(sum_error);
}

void output_progress(int iter, double error)
{
    std::string sentense;
    sentense += "\riter: ";
    sentense += std::to_string(iter);
    sentense += ", error: ";
    sentense += std::to_string(error);
    std::printf("%s", sentense.c_str());
}

void calculate()
{
    std::cout << "[SEQUENCE]: run multigrid" << std::endl;
    for(int iter = 1;;++iter)
    {
        // restriction
        int dim = std::pow(2, finest_power);
        gauss_seidel(0, dim);
        interp_restrict(0, dim);
        dim = dim/2;
        for(int level = 1; level <= maxlevel-1; ++level)
        {
            std::fill(u[level].begin(), u[level].end(), 0.0);
            gauss_seidel(level, dim);
            interp_restrict(level, dim);
            dim = dim/2;
        }
        gauss_seidel(maxlevel, dim); // now, dim = 2^(coarsest_power)
        
        // prolongation
        for(int level = maxlevel; level >= 1; --level)
        {
            interp_prolong(level, dim);
            dim *= 2;
        }

        double error = calculate_error();
        if(error < eps)
        {
            // finish!
            return;
        }

        output_progress(iter, error);
    }
}


int main()
{
    set_initial();
    calculate();
}
