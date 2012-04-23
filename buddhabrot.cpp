/*
 * Buddhabrot fractal generator
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


/* Usage: ./buddhabrot
 *
 * Set x_res, y_res, cpu_cores, max_orbit_size to suit your preferences.
 */
 
#include <iostream>
#include <complex>
#include <assert.h>
#include <cstdlib>
#include <cstdio>
#include <stdint.h>
#include "random.hpp"
#include <string.h>
#include <vector>

using namespace std;

/*
 * cpu_cores determines how many CPU cores the program will use.
 */
static const unsigned int cpu_cores = 3;

/*
 * x_res: horizontal resolution of any output images.
 * y_res: vertical resolution of any output images. 
 */
const unsigned int x_res = 4096;
const unsigned int y_res = 4096;

/*
 * max_orbit_size gives the maximum number of mandelbrot step iterations before
 * a value is considered to be a member of the mandelbrot set.
 */
const unsigned int max_orbit_size = 128;      
    
class BuddhabrotImager {

    /*
     * count stores how many times each pixel has been encountered.
     * +1 because floating point truncation is a bit dodgy
     */
    uint64_t           counts[x_res + 1][y_res + 1];
    /*
     * real/imag_min/max give the portion of the complex plane that are being
     * plotted.
     */
    float              real_min, real_max, imag_min, imag_max;

    /*
     * real_mult and imag_mult are values calculated from the above bounds. They
     * are used to speed up some floating point calculations. When the bounds
     * change, these are updated.
     */
    float              real_mult, imag_mult;

public:
    BuddhabrotImager() {}
    
    void setComplexBounds(float real_min, float real_max, float imag_min, float imag_max) {
        this->real_min = real_min;
        this->real_max = real_max;
        this->imag_min = imag_min;
        this->imag_max = imag_max;
        
        this->real_mult = x_res / (real_max - real_min);
        this->imag_mult = y_res / (imag_max - imag_min);
    }
    float  getRealMin() const {
        return real_min;
    }
    float  getRealMax() const {
        return real_max;
    }
    float  getImagMin() const {
        return imag_min;
    }
    float  getImagMax() const {
        return imag_min;
    }    
    
    /*
     * Set the entries of counts to zero.
     */
    void clearCounts(void) {    
        memset(this->counts, 0, sizeof(this->counts));
    }

    /*
     * Given a point in the complex plane, add one to the count for the corresponding pixel.
     * (i.e. slightly lightening that pixel). This should not be called for points that
     * are outside of the image area - see real_min, real_max, imag_min, imag_max.
     */
    void incrementCountFor(const complex<float>& c) {
        unsigned int x = (unsigned int)(real_mult * (c.real() - real_min));
        unsigned int y = (unsigned int)(imag_mult * (c.imag() - imag_min));
        
        assert(x <= x_res && y <= y_res);
        
        #pragma omp atomic    
        this->counts[x][y]++;
    }
    
    void convertPPMToPNG(const char* filename) const {
       char command[256];    
       FILE *fpipe;

       snprintf(command, 255, "convert %s %s.png", filename, filename);
       fpipe = popen(command, "r");
       pclose(fpipe);
       
       snprintf(command, 255, "rm %s", filename);
       fpipe = popen(command, "r");
       pclose(fpipe);   
    }

    /*
     * Counts must be scaled before this function is called.
     */
    void producePNG(unsigned int file_id) const {
        char  filename[256];
        char *write_buffer = (char *)malloc(1000000);
        FILE *f;
        
        snprintf(filename, sizeof(filename), "images/buddhabrot-%u.ppm", file_id);
        f = fopen(filename, "w");
        
        /*
         * Enabling buffering should speed up writing the file.
         */
        setbuf(f, write_buffer);
               
        fprintf(f, "P3\n%u\n%u\n255\n", x_res, y_res);
        for (unsigned int y = 0; y < y_res; y++) {
            for (unsigned int x = 0; x < x_res; x++) {
                unsigned int out = counts[x][y];
                fprintf(f, "%u %u %u\n", out, out, out);
            }
        }
        fclose(f);
        
        convertPPMToPNG(filename);
        free(write_buffer);
    }

    bool isInRange(const complex<float>& z) const {
        return (z.real() >= real_min &&
                z.real() <= real_max &&
                z.imag() >= imag_min &&
                z.imag() <= imag_max);
    }
    
    /*
     * Find an area with a good range of brightness values and adjust the bounds to zoom into it.
     */
    void zoomIntoInterestingArea(void)
    {
        unsigned int step_size = 128;
        unsigned int zoom_size = 512;
        unsigned int start_x, start_y;
        
        unsigned int histogram[256];

        unsigned int best_start_x = 0, best_start_y = 0;
        float        best_entropy = 0.0f;
        
        /*
         * Try lots of possible regions. For each region, count the number of pixel brightnesses
         */
        for (start_x = 0; start_x <= x_res - zoom_size; start_x += step_size) {
            for (start_y = 0; start_y <= y_res - zoom_size; start_y += step_size) {
                                
                for (unsigned int i = 0; i < 256; i++) {
                    histogram[i] = 0;
                }
                
                for (unsigned int x = start_x; x < start_x + zoom_size; x++) {
                    for (unsigned int y = start_y; y < start_y + zoom_size; y++) {
                        assert(counts[x][y] < 256);
                        histogram[counts[x][y]]++;
                    }
                }
                
                float entropy = 0.0f;
                for (unsigned int i = 0; i < 256; i++) {
                    if (histogram[i] > 0) {
                        float p = histogram[i] / (float)(zoom_size * zoom_size);
                        entropy -= p * log2f(p);
                    }
                }
                
                if (entropy > best_entropy) {
                    best_entropy = entropy;
                    best_start_x = start_x;
                    best_start_y = start_y;
                }
            }
        }
        
        std::cout << "Best entropy: " << best_entropy << "\n";
        
        float new_r_min = real_min + (real_max - real_min) * 
                                                     (best_start_x / (float)x_res);
        float new_r_max = real_min + (real_max - real_min) *
                                       ((best_start_x + zoom_size) / (float)x_res);
        float new_i_min = imag_min + (imag_max - imag_min) * 
                                                     (best_start_y / (float)y_res);
        float new_i_max = imag_min + (imag_max - imag_min) *
                                       ((best_start_y + zoom_size) / (float)y_res);
            
        /*
         * Use setComplexBounds to update real_mult and imag_mult at the same time.
         */
        this->setComplexBounds(new_r_min, new_r_max, new_i_min, new_i_max);
    }
    
    void scaleCounts(uint64_t  max_value)    
    {
        uint64_t total = 0;
        
        /*
         * Calculate the average pixel value by averaging all the rows and then averaging that.
         */
        for (unsigned int x = 0; x < x_res; x++) {
            uint64_t x_total = 0;
            for (unsigned int y = 0; y < y_res; y++) {
                x_total += counts[x][y];
            }
            total += x_total / y_res;
        }
        
        uint64_t average = total / x_res;
        
        /*
         * Avoid division by zero.
         */
        if (average == 0) {
            average = 1;
        }
        
        for (unsigned int y = 0; y < y_res; y++) {
            for (unsigned int x = 0; x < x_res; x++) {
                /*
                 * Somewhat arbitrary scaling that seems to work ok.
                 */
                unsigned int out = (counts[x][y] * max_value) / (5 * average);
                if (out > max_value) {
                    out = max_value;
                }
                
                counts[x][y] = out;
            }
        }
    }
};

/*
 * Returns |z|^2.
 */
static inline float
magnitude_squared (const complex<float>& z)
{
    return (z.imag() * z.imag()) + (z.real() * z.real());
}


/*
 * DistributionGenerator uses the Metropolis-Hastings Algorithm to generate
 * random samples x uniformly from the subset of D(2) (the open complex disc of
 * radius 2) with orbit(x) intersect {the region of the image being examined} non-empty.
 */
struct DistributionGenerator
{
private:
    complex<float>            previous_sample;    
    uint64_t                  unrepeated_sample_count;
    bool                      any_found;
    float                     sigma;
    RNG                       rng;
    bool                      last_sample_was_new;
    
    /* 
     * Space to store the orbit for the previous sample and a new proposed sample.
     */
    vector< complex<float> >  orbit[2];   
    unsigned int              orbit_id;
    
    /*
     * Uniform distribution is faster but probably worse.
     */
#ifdef NORMAL_DISTRIBUTION
    RandNormal2D              random_offset;
#else
    RandUniform2D             random_offset;
#endif

    BuddhabrotImager          &imager;

public:
    DistributionGenerator(BuddhabrotImager  &_imager) : imager(_imager)
    {
        unrepeated_sample_count = 0; 
        sigma                   = 0.05f;
        previous_sample         = complex<float>(2.0f, 2.0f);
        orbit_id                = 0;
        any_found               = false;
        
        orbit[0].reserve(max_orbit_size);
        orbit[1].reserve(max_orbit_size);        
    }
    
    const vector< complex<float> >& getOrbit() 
    {
        return orbit[orbit_id];
    }
    
    bool isOrbitInRange(complex<float> c)
    {
        float          z_size      = magnitude_squared(c);
        complex<float> z           = c;       

        orbit[1 - orbit_id].clear();
                
        for (unsigned int i = 0; i < max_orbit_size && z_size <= 4.0f; i++) {
            if (imager.isInRange(z)) {                
                orbit[1 - orbit_id].push_back(z);
            }
            
            z *= z;
            z += c;
            z_size = magnitude_squared(z);
        }

        return (orbit[1 - orbit_id].size() > 0) && (z_size > 4);
    }
    
    bool findStartingValue(uint64_t  max_iterations)
    {       
        bool  any_samples_found = false;

        for (uint64_t  i = 0; i < max_iterations && !any_samples_found; i++) {
            complex<float> sample = complex<float>(rng() * 4.0f - 2.0f,
                                                   rng() * 4.0f - 2.0f);
            any_samples_found = this->isOrbitInRange(sample);
        }
        return (any_samples_found);
    }
    
    bool wasLastSampleNew() const
    {
        return this->last_sample_was_new;
    }
    
    void iterate(void)
    {
        complex<float> proposed_new_sample;
        const float    chance_random  = 0.1;
        
        /*
         * Generate a proposed new sample value.
         * - With chance 'chance_random' we choose a random value in the complex disc of radius 2.
         * - Otherwise we choose a value close to the previous sample. A 2-dimensional normal
         *   distribution sample is added to the previous sample. The standard deviation of this
         *   normal distribution is determined by sigma.        
         */
        if (this->rng() < chance_random) {
            do {
                proposed_new_sample = complex<float>(rng() * 4.0f - 2.0f,
                                                     rng() * 4.0f - 2.0f);
            } while (magnitude_squared(proposed_new_sample) >= 4.0f);
            
        } else {
            proposed_new_sample = previous_sample + random_offset() * sigma;           
        }
     
        if (this->isOrbitInRange(proposed_new_sample)) {
            last_sample_was_new = true;
            previous_sample =  proposed_new_sample;

            orbit_id = 1 - orbit_id;          
        } else {
            last_sample_was_new = false;
        }
    }
};

void createBuddhabrotImages(BuddhabrotImager &bi)
{

    /*
     * Starting bounds.
     */
    bi.setComplexBounds(-2.0f, 2.0f, -2.0f, 2.0f);
  
    /*
     * We don't set a maximum zoom - we'll quit when the assert statement a bit
     * later fails.
     */
    for (unsigned int zoom = 1;; zoom++) {    
        std::cout << "Zoom level " << zoom
                  << " real: "  << bi.getRealMin() <<  " to " << bi.getRealMax()
                  << ", imag: " << bi.getImagMin() << "i to " << bi.getImagMax() << "i \n";
        bi.clearCounts();       
        
        #pragma omp parallel for
        for (unsigned int core = 0; core < cpu_cores; core++) {       
            DistributionGenerator   dg(bi);

            if (!dg.findStartingValue(500000)) {
                #pragma omp critical
                {
                    printf("Couldn't find even 1 value in the zoom window.\n");
                    exit(1);
                }
            }

            unsigned int rep_max           = 200000000UL;
            unsigned int rep_percent_jump  = rep_max / 1000;
            unsigned int rep_print_percent = rep_percent_jump;        
            unsigned int start             = time(NULL);
            uint64_t     genuine_samples   = 0;
            
            for (unsigned int rep = 0; rep < rep_max; rep++) {
                dg.iterate();
                
                if (dg.wasLastSampleNew()) {                   
                    genuine_samples++;
                }
                
                /*
                 * Output percentage completion statistics on the first core.
                 */
                if (core == 0 && rep == rep_print_percent) {
                    float percent = (100ULL * rep) / float(rep_max);
                    float eta     = (float(time(NULL) - start) * (rep_max - rep)) / rep;
                    
                    std::cerr << percent << "% (eta " << (eta / 60) << " mins)" 
                              << "                  " << std::flush << "\r";
                    
                    rep_print_percent += rep_percent_jump;
                }
                
                /*
                 * Only points that are in range are included in the orbit.
                 */
                const vector< complex<float> > &orbit = dg.getOrbit();
                for (unsigned int i = 0; i < orbit.size(); i++) {                   
                    bi.incrementCountFor(orbit[i]);
                }
            }
            
            #pragma omp critical
            {
                std::cerr << (100.0f * float(genuine_samples) / rep_max) 
                          << "% good samples                   \n";            
                std::cerr << "thread " << core << " done in " 
                          << (time(NULL) - start) << " seconds       \n";                          
            }
        }

        bi.scaleCounts(255);
        bi.producePNG(zoom);
        bi.zoomIntoInterestingArea();
    }
}

int main(int argc, char* argv[])
{
    /*
     * Create on heap as stack not large enough.
     */
    BuddhabrotImager *bi_ptr = new BuddhabrotImager();     
    createBuddhabrotImages(*bi_ptr);    
    delete bi_ptr;
}
