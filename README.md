# SimpleBlueNoise
A small, header-only blue noise texture generator in C99.

# Usage
To use the library, simply include the header file. You must **#define** the implementation in exactly _one_ .c/.cpp file.

```c
#define BLUE_NOISE_GENERATOR_IMPLEMENTATION
#include "blue_noise_generator.h"
```

The library includes two different algorithms to produce blue noise textures. Both functions expect a pointer to a contiguous array of **unsigned int**, the size of which should be equal to **width** × **height**. 

```c
int blue_noise_generator_create_void_and_cluster(unsigned int* buffer, int width, int height);
int blue_noise_generator_create_forced_random(unsigned int* buffer, int width, int height);
```

Each function returns 1 on success and 0 on failure. If successful, every value in the buffer will be given a unique rank starting from 0 counting up to **width** × **height** - 1. This should make it easy to normalise the values to fit your desired texture format.

The library uses its own psuedorandom number generator, and includes a function to set the current seed:

```c
void blue_noise_generator_set_seed(unsigned int seed);
```

# What is Blue Noise?
Blue noise is noise that contains mostly high-frequency components, meaning that the value difference between any two adjacent pixels is generally greater than that of white noise. This gives it a finer, more granular quality, as opposed to white noise which can appear more coarse.

| White Noise | Blue Noise |
|-|-|
|![white_noise_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/5eeb3dca-8d83-4853-b694-455c62ae9922)|![void_and_cluster_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/f30335d9-010b-4ade-b3df-122f01afceb0)|

Blue noise textures are useful whenever you need evenly distributed random samples. Some common use cases include image dithering and 3D rendering.

| White Noise Dithering | Blue Noise Dithering |
|-|-|
|![dither_white_noise](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/d66d5bd7-bc75-4666-8257-a0d0181d4c1e)|![dither_blue_noise](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/0e7f585e-ba96-4306-8ed5-7232e68b228b)|

# Implementation
### Void-and-Cluster Algorithm
This library uses a variant of Robert Ulichney's Void-and-Cluster method for dither array generation[^1]. It makes use of a look-up table to simplify energy computations. The core of the algorithm involves evaluating the following Gaussian energy function:

$$\Huge{f(x) = e^{-{{r}^2}/{2\sigma^2}}}$$

Where $\large{r}$ is the wrap-around (toroidal) distance between pixels. The parameter $\large{\sigma}$ can be set via the following function:

```c
void blue_noise_generator_set_void_and_cluster_sigma(double sigma);
```

When the binary pattern is first created, a number of pixels are randomly selected as minority pixels to set the initial state. The number of minority pixels can be controlled by a threshold parameter, which should range from 0.0 to 1.0:

```c
void blue_noise_generator_set_void_and_cluster_threshold(double threshold);
```

### Forced Random Algorithm
This algorithm is taken from the paper _"Forced random dithering: improved threshold matrices for ordered dithering"_[^2]. It is less thorough than the void-and-cluster method but is faster as a result. It uses the following energy function:

$$\Huge{f(x) = e^{-({r}/{d})^s}}$$

The parameters $\large{d}$ and $\large{s}$ control the deviation and steepness of the normal curve, respectively. They can be set via the following functions:

```c
void blue_noise_generator_set_forced_random_deviation(double deviation);
void blue_noise_generator_set_forced_random_steepness(double steepness);
```

Every iteration, a random number of candidate pixels are selected and compared to determine the next sample point. The number of candidate pixels can be controlled by a threshold parameter, which should range from 0.0 to 1.0:

```c
void blue_noise_generator_set_forced_random_threshold(double threshold);
```

# Example Output

Below is a comparison between the output of each algorithm and their corresponding [spectral densities](https://en.wikipedia.org/wiki/Spectral_density). The energy parameters have been set so that the normal curve given by each energy function is approximately the same.

|Void-and-Cluster, σ = 1.9, threshold = 0.1|
|-|
|![void_and_cluster_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/f30335d9-010b-4ade-b3df-122f01afceb0) ![void_and_cluster_256_fft](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/735a57b8-3c01-406b-82e2-3f9e3e3590f6)|

|Forced Random, s = 2.0, d = 2.7, threshold = 0.5 |
|-|
|![forced_random_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/763280c9-6ea2-4042-96e4-8513e915818c) ![forced_random_256_fft](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/bd6cc5c6-92fc-4309-9f79-dc89e56b27f3)|

# Performance

The table below shows the run time of each algorithm for various texture resolutions, taken as an average of 10 iterations. We can see that Forced Random is about twice as fast as Void-and-Cluster.

Platform: _Intel(R) Core(TM) i7-8850H CPU @ 2.60GHz_.

| Resolution | Void-and-Cluster | Forced Random |
|:-:         | :-:              | :-:           |
| 16 x 16    | 2ms      | <1ms    |
| 32 x 32    | 15ms     | 8ms     |
| 64 x 64    | 231ms    | 128ms   |
| 128 x 128  | 4s       | 2s      |
| 256 x 256  | 62s      | 37s     |

[^1]: R. Ulichney, _"The void-and-cluster method
for dither array generation"_ (1993).

[^2]: W. Purgathofer, R.F. Tobler & M. Geiler, _"Forced random dithering: improved threshold matrices for ordered dithering"_ (1994).
