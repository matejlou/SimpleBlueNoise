#ifndef BLUE_NOISE_GENERATOR_H
#define BLUE_NOISE_GENERATOR_H

#ifdef __cplusplus
extern "C" {
#endif

#ifndef BLUE_NOISE_GENERATOR_LINKAGE
    #ifdef BLUE_NOISE_GENERATOR_STATIC
        #define BLUE_NOISE_GENERATOR_LINKAGE static
    #else
        #define BLUE_NOISE_GENERATOR_LINKAGE extern
    #endif
#endif

/* Create blue noise using the void-and-cluster algorithm. Values in the buffer start from 0 and end at (width * height - 1). */
BLUE_NOISE_GENERATOR_LINKAGE int blue_noise_generator_create_void_and_cluster(unsigned int* buffer, int width, int height);

/* Create blue noise using the forced random algorithm. Values in the buffer start from 0 and end at (width * height - 1). */
BLUE_NOISE_GENERATOR_LINKAGE int blue_noise_generator_create_forced_random(unsigned int* buffer, int width, int height);

/* Set the sigma term of the energy computation for the void-and-cluster algorithm. */
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_void_and_cluster_sigma(double sigma);

/* Set the initial number of random samples as a proportion of total samples from 0.0 to 1.0. */
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_void_and_cluster_threshold(double threshold);

/* Set the steepness of the energy computation for the forced random algorithm. */
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_steepness(double steepness);

/* Set the deviation of the energy computation for the forced random algorithm. */
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_deviation(double deviation);

/* Set the number of samples tested during each iteration as a proportion of available samples from 0.0 to 1.0.*/
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_threshold(double threshold);

/* Set the seed for the random number generator. */
BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_seed(unsigned int seed);

#ifdef __cplusplus
}
#endif

#endif // !BLUE_NOISE_GENERATOR_H

#ifdef BLUE_NOISE_GENERATOR_IMPLEMENTATION

    /* Implementation starts here. */

#include <math.h> /* exp, INFINITY */
#include <string.h> /* memset, memcpy */

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

#if defined(BLUE_NOISE_GENERATOR_INTERNAL_MALLOC) != defined(BLUE_NOISE_GENERATOR_INTERNAL_FREE)
#error "Must define both or none of BLUE_NOISE_GENERATOR_INTERNAL_MALLOC and BLUE_NOISE_GENERATOR_INTERNAL_FREE."
#endif

#ifndef BLUE_NOISE_GENERATOR_INTERNAL_MALLOC
    #include <malloc.h> /* malloc, free */
    #define BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(size) malloc(size)
    #define BLUE_NOISE_GENERATOR_INTERNAL_FREE(ptr) free(ptr)
#endif

static const unsigned int BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_MAX = 0xffff;
static unsigned int BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_STATE = 1;

    /* VOID-AND-CLUSTER ALGORITHM */

struct BlueNoiseGeneratorInternalVoidAndClusterData
{
    int width, height;
    unsigned int* buffer;
    double* energy_table;
    double* energy_kernel;
    struct { int x, y; double energy; } highest, lowest;
};

enum BlueNoiseGeneratorInternalEnum
{   
    BLUE_NOISE_GENERATOR_INTERNAL_ADD = 1,
    BLUE_NOISE_GENERATOR_INTERNAL_SUB = -1,
};

static double BLUE_NOISE_GENERATOR_INTERNAL_VAC_SIGMA = 1.9;
static double BLUE_NOISE_GENERATOR_INTERNAL_VAC_THRESHOLD = 0.1;

static void blue_noise_generator_internal_void_and_cluster_update_energy(struct BlueNoiseGeneratorInternalVoidAndClusterData* data, const int o_x, const int o_y, const int operation, const int invert)
{
    // Reset highest/lowest energy
    data->highest.energy = -INFINITY;
    data->lowest.energy = +INFINITY;

    for (int e_y = 0; e_y < data->height; e_y++)
    {
        // Get the toroidal y distance from the origin
        int dy = abs(o_y - e_y);
        if (dy > data->height / 2) dy = data->height - dy;

        const double energy_y = data->energy_kernel[dy];
        
        for(int e_x = 0; e_x < data->width; e_x++)
        {
            // Get the toroidal x distance from the origin
            int dx = abs(o_x - e_x);
            if (dx > data->width / 2) dx = data->width - dx;

            // Increment or decrement by the value in the energy mask
            const double energy_x = data->energy_kernel[dx];
            data->energy_table[e_x + e_y * data->width] += energy_x * energy_y * (double)operation;

            // Update highest and lowest values
            if (invert)
            {
                // If we're in phase 3 (invert meaning of voids and clusters) then we only care about the highest energy
                if (data->energy_table[e_x + e_y * data->width] > data->highest.energy &&
                    data->buffer[e_x + e_y * data->width] < 1)
                {
                    data->highest.energy = data->energy_table[e_x + e_y * data->width];
                    data->highest.x = e_x;
                    data->highest.y = e_y;
                }
            }
            else
            {
                if (data->energy_table[e_x + e_y * data->width] > data->highest.energy &&
                    data->buffer[e_x + e_y * data->width] > 0)
                {
                    data->highest.energy = data->energy_table[e_x + e_y * data->width];
                    data->highest.x = e_x;
                    data->highest.y = e_y;
                }
                if (data->energy_table[e_x + e_y * data->width] < data->lowest.energy &&
                    data->buffer[e_x + e_y * data->width] < 1)
                {
                    data->lowest.energy = data->energy_table[e_x + e_y * data->width];
                    data->lowest.x = e_x;
                    data->lowest.y = e_y;
                }
            }
        }
    }
}

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_void_and_cluster_sigma(const double sigma)
{
    BLUE_NOISE_GENERATOR_INTERNAL_VAC_SIGMA = sigma;
}

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_void_and_cluster_threshold(const double threshold)
{
    BLUE_NOISE_GENERATOR_INTERNAL_VAC_THRESHOLD = threshold;
}

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_seed(const unsigned int seed)
{
    BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_STATE = seed;
}

BLUE_NOISE_GENERATOR_LINKAGE unsigned int blue_noise_generator_internal_random()
{
    BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_STATE = 214013u * BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_STATE + 2531011u;
    return (BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_STATE >> 16) & BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_MAX;
}

BLUE_NOISE_GENERATOR_LINKAGE int blue_noise_generator_create_void_and_cluster(unsigned int* buffer, const int width, const int height)
{
    // Allocate the energy table and energy mask
    const int half_width = width / 2 + 1;
    const int half_height = height / 2 + 1;
    const int kernel_size = (width > height) ? half_width : half_height; // The space is toroidal, so we only need half the size
    double* energy_table = (double*)BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(sizeof(*energy_table) * width * height);
    double* energy_kernel = (double*)BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(sizeof(*energy_kernel) * kernel_size); 

    if (energy_table == NULL || energy_kernel == NULL) 
    {
        if (energy_table != NULL) BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table);
        if (energy_kernel != NULL) BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_kernel);
        return 0;
    }

    memset(energy_table, 0, sizeof(*energy_table) * width * height);

    // Initialise data struct
    struct BlueNoiseGeneratorInternalVoidAndClusterData data;
    data.width          = width;
    data.height         = height;
    data.buffer         = buffer;
    data.energy_table   = energy_table;
    data.energy_kernel  = energy_kernel;
    data.highest.energy = -INFINITY;
    data.lowest.energy  = +INFINITY;

    // Compute the energy kernel. Because our energy filter is separable, we only need to precompute the energy values along one axis! :)
    for (int i = 0; i < kernel_size; i++)
    {
        const double distance = (double)(i * i);
        const double energy = exp(-distance / (2.0 * BLUE_NOISE_GENERATOR_INTERNAL_VAC_SIGMA * BLUE_NOISE_GENERATOR_INTERNAL_VAC_SIGMA));
        energy_kernel[i] = energy;
    }

    // Create white noise input points
    unsigned int sample_count = 0;

    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            if ((double)blue_noise_generator_internal_random() / (double)BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_MAX < BLUE_NOISE_GENERATOR_INTERNAL_VAC_THRESHOLD)
            {
                blue_noise_generator_internal_void_and_cluster_update_energy(&data, x, y, BLUE_NOISE_GENERATOR_INTERNAL_ADD, 0);
                buffer[x + y * width] = 1;
                sample_count += 1;
            }
            else
            {
                buffer[x + y * width] = 0;
            }
        }
    }

    // Process input points to create the Prototype Binary Pattern
    while (1)
    {
        buffer[data.highest.x + data.highest.y * width] = 0;
        blue_noise_generator_internal_void_and_cluster_update_energy(&data, data.highest.x, data.highest.y, BLUE_NOISE_GENERATOR_INTERNAL_SUB, 0);

        const int minimum_x = data.lowest.x;
        const int minimum_y = data.lowest.y;

        buffer[data.lowest.x + data.lowest.y * width] = 1;
        blue_noise_generator_internal_void_and_cluster_update_energy(&data, data.lowest.x, data.lowest.y, BLUE_NOISE_GENERATOR_INTERNAL_ADD, 0);

        if (minimum_x == data.highest.x && minimum_y == data.highest.y)
            break;
    }

    /* PHASE 1 */

    // Make the sample points progressive
    double* energy_table_temp = (double*)malloc(sizeof(*energy_table_temp) * width * height);

    if (energy_table_temp == NULL) 
    {
        BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table);
        BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_kernel);
        return 0;
    }

    memcpy(energy_table_temp, energy_table, sizeof(*energy_table) * width * height);
    data.energy_table = energy_table_temp; // Change pointer

    unsigned int sample_count_temp = sample_count;

    while (sample_count_temp > 0)
    {
        buffer[data.highest.x + data.highest.y * width] = sample_count_temp--;
        blue_noise_generator_internal_void_and_cluster_update_energy(&data, data.highest.x, data.highest.y, BLUE_NOISE_GENERATOR_INTERNAL_SUB, 0);
    }

    BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table_temp);
    data.energy_table = energy_table; // Restore pointer to original table

    /* PHASE 2 */

    // Add samples until half the pattern is filled
    while (sample_count < (unsigned int)(width * height / 2))
    {
        buffer[data.lowest.x + data.lowest.y * width] = ++sample_count;
        blue_noise_generator_internal_void_and_cluster_update_energy(&data, data.lowest.x, data.lowest.y, BLUE_NOISE_GENERATOR_INTERNAL_ADD, 0);
    }
    
    /* PHASE 3 */

    // Recompute the energy table such that voids emit energy
    memset(energy_table, 0, sizeof(*energy_table) * width * height);

    for (int y = 0; y < height; y++)  
    {
        for (int x = 0; x < width; x++) 
        {
            if (buffer[x + y * width] > 0)
                continue;

            blue_noise_generator_internal_void_and_cluster_update_energy(&data, x, y, BLUE_NOISE_GENERATOR_INTERNAL_ADD, 1);
        }
    }

    // Continue to add samples until the whole pattern is filled
    while (sample_count < (unsigned int)(width * height))
    {
        buffer[data.highest.x + data.highest.y * width] = ++sample_count;
        blue_noise_generator_internal_void_and_cluster_update_energy(&data, data.highest.x, data.highest.y, BLUE_NOISE_GENERATOR_INTERNAL_SUB, 1);
    }

    // Subtract 1 from buffer values so that the starting value is 0
    for (int i = 0; i < width * height; i++)
    {
        buffer[i] -= 1;
    }

    BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table);
    BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_kernel);

    return 1;
}

    /* FORCED RANDOM ALGORITHM */

struct BlueNoiseGeneratorInternalForcedRandomCandidate
{
    int x, y;
};

static double BLUE_NOISE_GENERATOR_INTERNAL_FR_STEEPNESS = 2.0;
static double BLUE_NOISE_GENERATOR_INTERNAL_FR_DEVIATION = 2.7;
static double BLUE_NOISE_GENERATOR_INTERNAL_FR_THRESHOLD = 0.5;

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_steepness(const double steepness)
{
    BLUE_NOISE_GENERATOR_INTERNAL_FR_STEEPNESS = steepness;
}

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_deviation(const double deviation)
{
    BLUE_NOISE_GENERATOR_INTERNAL_FR_DEVIATION = deviation;
}

BLUE_NOISE_GENERATOR_LINKAGE void blue_noise_generator_set_forced_random_threshold(const double threshold)
{
    BLUE_NOISE_GENERATOR_INTERNAL_FR_THRESHOLD = threshold;
}

BLUE_NOISE_GENERATOR_LINKAGE int blue_noise_generator_create_forced_random(unsigned int* buffer, const int width, const int height)
{
    // Allocate the energy table and energy mask
    const int half_width = width / 2 + 1;
    const int half_height = height / 2 + 1;
    double* energy_table = (double*)BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(sizeof(*energy_table) * width * height);
    double* energy_kernel = (double*)BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(sizeof(*energy_kernel) * half_width * half_height);
    struct BlueNoiseGeneratorInternalForcedRandomCandidate* candidates = (struct BlueNoiseGeneratorInternalForcedRandomCandidate*)BLUE_NOISE_GENERATOR_INTERNAL_MALLOC(sizeof(*candidates) * width * height);

    if (energy_table == NULL || energy_kernel == NULL || candidates == NULL) 
    {
        if (energy_table != NULL) BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table);
        if (energy_kernel != NULL) BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_kernel);
        if (candidates != NULL) BLUE_NOISE_GENERATOR_INTERNAL_FREE(candidates);
        return 0;
    }

    memset(buffer, 0, sizeof(*buffer) * width * height);
    memset(energy_table, 0, sizeof(*energy_table) * width * height);

    // Fill the candidate array
    for (int y = 0; y < height; y++)
    {
        for (int x = 0; x < width; x++)
        {
            candidates[x + y * width].x = x;
            candidates[x + y * width].y = y;
        }
    }

    // Compute the energy kernel. This energy function is not separable (?) so we must compute the full 2D quadrant
    for (int y = 0; y < half_height; y++)
    {
        for (int x = 0; x < half_width; x++)
        {
            const double distance = sqrt((double)(x * x + y * y));
            const double energy = exp(-pow(distance / BLUE_NOISE_GENERATOR_INTERNAL_FR_DEVIATION, BLUE_NOISE_GENERATOR_INTERNAL_FR_STEEPNESS));
            energy_kernel[x + y * half_width] = energy;
        }
    }

    // Progressively add sample points until the pattern is filled
    unsigned int sample_count = 0;
    struct { double energy; int candidate; } minimum = { 0.0, 0 };
    
    while (sample_count < (unsigned int)(width * height))
    {
        minimum.energy = INFINITY;

        // Select n candidates at random from the available pool of candidates
        int number_of_available_candidates = width * height - sample_count;
        int n = (int)((double)(number_of_available_candidates) * BLUE_NOISE_GENERATOR_INTERNAL_FR_THRESHOLD);
        
        if (n == 0) 
        {
            n = 1;
        }

        // Shuffle n candidates, keeping track of the lowest energy
        for (int i = 0; i < n; i++)
        {
            const int index = (int)((double)blue_noise_generator_internal_random() / ((double)BLUE_NOISE_GENERATOR_INTERNAL_RANDOM_MAX + 1.0) * number_of_available_candidates);
            const int x = candidates[index].x;
            const int y = candidates[index].y;
            const double current_energy = energy_table[x + y * width];

            // Put the selected candidate at the back of the array
            const struct BlueNoiseGeneratorInternalForcedRandomCandidate temp = candidates[index];
            candidates[index] = candidates[number_of_available_candidates - 1];
            candidates[number_of_available_candidates - 1] = temp;

            if (current_energy < minimum.energy)
            {
                minimum.energy = current_energy;
                minimum.candidate = number_of_available_candidates - 1;
            }

            number_of_available_candidates--;
        }

        // Add the next sample to the minimum position and accummulate the energy mask
        const int minimum_x = candidates[minimum.candidate].x;
        const int minimum_y = candidates[minimum.candidate].y;
        buffer[minimum_x + minimum_y * width] = ++sample_count;

        // Put the minimum candidate at the back of the array
        const int last_candidate = width * height - sample_count;
        const struct BlueNoiseGeneratorInternalForcedRandomCandidate temp = candidates[minimum.candidate];
        candidates[minimum.candidate] = candidates[last_candidate];
        candidates[last_candidate] = temp;
        
        for (int y = 0; y < height; y++)
        {
            // Get the toroidal y distance from the origin
            int dy = abs(minimum_y - y);
            if (dy > height / 2) dy = height - dy;

            for(int x = 0; x < width; x++)
            {
                // Get the toroidal x distance from the origin
                int dx = abs(minimum_x - x);
                if (dx > width / 2) dx = width - dx;

                // Increment or decrement by the value in the energy mask
                const double energy = energy_kernel[dx + dy * half_width];
                energy_table[x + y * width] += energy;
            }
        }
    }

    // Subtract 1 from buffer values so that the starting value is 0;
    for (int i = 0; i < width * height; i++)
    {
        buffer[i] -= 1;
    }

    BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_table);
    BLUE_NOISE_GENERATOR_INTERNAL_FREE(energy_kernel);
    BLUE_NOISE_GENERATOR_INTERNAL_FREE(candidates);

    return 1;
}

#endif

/*
    MIT License

    Copyright (c) 2023 matejlou

    Permission is hereby granted, free of charge, to any person obtaining a copy
    of this software and associated documentation files (the "Software"), to deal
    in the Software without restriction, including without limitation the rights
    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
    copies of the Software, and to permit persons to whom the Software is
    furnished to do so, subject to the following conditions:

    The above copyright notice and this permission notice shall be included in all
    copies or substantial portions of the Software.

    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
    SOFTWARE.
*/
