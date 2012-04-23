/*
 * Random number test module, part of buddhabrot fractal generator
 * Copyright (C) 2010-2012 Ryan Lothian
 *
 * This program is free software; you can redistribute it and/or modify it
 * under the terms of the GNU General Public License as published by the Free
 * Software Foundation; either version 2 of the License, or (at your option) 
 * any later version.
 *
 * This program is distributed in the hope that it will be useful, but WITHOUT
 * ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
 * FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for 
 * more details.
 *
 * You should have received a copy of the GNU General Public License along with
 * this program; if not, write to the Free Software Foundation, Inc., 
 * 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA
 */
 
#include "random.hpp"
#include <iostream>
#include <cstdio>
#include <cstdlib>
#include <assert.h>

#undef NDEBUG

int main(int argc, char* argv[])
{
    RNG rng;

    /*
     * Uniform RNG.
     *
     * Save the outputted values to a CSV. The user should check them by hand.
     * Suggested checks: average ~= 0.5, stdev
     *                   histogram with splits at 0.0, 0.1, 0.2, ..., 1.0
     */
    FILE *csv_out = fopen("uniform.csv", "w");
    for (unsigned int i = 0; i < 1000; i++) {
        float r = rng();
        
        assert(0.0f <= r && r <= 1.0f);
        
        fprintf(csv_out, "%g\n", r);
    }
    fclose(csv_out);
    
    
    /*
     * Normal RNG.
     
     * Save the outputted values to a CSV. The user should check them by hand.
     * Suggested checks: average ~= 0.0, std. dev ~= 1.0
     */     
    RandNormal2D normal;
    
    csv_out = fopen("normal.csv", "w");
    for (unsigned int i = 0; i < 1000; i++) {
        std::complex<float> r = normal();
        fprintf(csv_out, "%g,%g\n", r.real(), r.imag());
    }
    fclose(csv_out);
    
    std::cout << "The user should now check normal.csv's randomness.\n";

    return 0;
}
