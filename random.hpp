/*
 * Random number generator module, part of buddhabrot fractal generator
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
 
#include <cstdlib>
#include <stdint.h>
#include <complex>
#include <omp.h>

static const float TWO_PI = 2.0f * 3.14159265358979f;

/*
 * RNG generates random floats in the closed interval [0, 1] using an LCG. It should be faster than
 * rand() (as there is no library call) and produce at least as good randomness (64-bit LCG).
 */
struct RNG {
    uint64_t state;
    
    RNG() {
        /*
         * Seed with rand() and omp_get_thread_num().
         */
        #pragma omp critical
        {
            state = ((uint64_t)rand() << 32) ^ (rand() + omp_get_thread_num());
        }
    }
    
    /*
     * Generate a floating-point value in [0, 1].
     */
    float operator () () {
        state *= 6364136223846793005ULL;
        state += 1442695040888963407ULL;
        return float(state) / float((uint64_t)(-1LL));
    }
};


/*
 * RandNormal2D generates complex numbers whose real and imaginary components
 * are both drawn from a normal distribution with mean 0 and s.d. 1.
 */
struct RandNormal2D {
    RNG rng;
    
    RandNormal2D() {}
    
    std::complex<float> operator () () {
        float u1        = rng();
        float u2        = rng();
        float r         = sqrtf(-2.0f * logf(u1));
        float theta     = TWO_PI * u2;
        float cos_theta = cosf(theta);
        float sin_theta = sinf(theta);

        return std::complex<float>(r * cos_theta, r * sin_theta);
    }
};

/*
 * RandUniform2D generates complex numbers whose real and imaginary components
 * are both drawn from a uniform distribution with mean 0 and range 1.
 */
    RNG rng;
    
struct RandUniform2D {
    RandUniform2D() {}
    
    std::complex<float> operator () () {
        return std::complex<float>((rng() * 2.0f) - 1.0f,
                                   (rng() * 2.0f) - 1.0f);
    }
};

