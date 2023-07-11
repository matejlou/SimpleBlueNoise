# SimpleBlueNoise
A small, header-only blue noise texture generator in C99.

# Usage
To use the library, simply **#include** the header file. You must **#define** the implementation in exactly _one_ .c/.cpp file.

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
|![white_noise_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/eda47b5f-cd89-4195-88a4-9f108eaabe44)|![void_and_cluster_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/c85d399e-2ab6-4c18-92cb-b3495847ab01)|

Blue noise textures are useful whenever you need evenly distributed random samples. Some common use cases include image dithering and 3D rendering.

| White Noise Dithering | Blue Noise Dithering |
|-|-|
|![dither_white_noise](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/c771a2aa-c511-4937-b37c-b4069a2bb18e)|![dither_blue_noise](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/6356884a-5a65-49d1-bfa8-bcebe10c4b7d)|

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
|![void_and_cluster_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/da03258d-b112-4b26-8e5c-a383d61b4aee) ![void_and_cluster_256_fft](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/19586558-daca-4b83-8bec-c36f8a05cebf)|

|Forced Random, s = 2.0, d = 2.7, threshold = 0.5 |
|-|
|![forced_random_256](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/b9de18d3-91d3-41ec-a751-b4c1e9380d01) ![forced_random_256_fft](https://github.com/matejlou/SimpleBlueNoise/assets/120740455/961b5c83-14f8-489e-bafe-de54fec99557)|

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

# Additional Reading

[Christoph Peters, _"Free Blue Noise Textures."_](http://momentsingraphics.de/BlueNoise.html)

[Alan Wolfe, _"What the Heck is Blue Noise?"_](https://blog.demofox.org/2018/01/30/what-the-heck-is-blue-noise/)

[Alan Wolfe, _"Using Blue Noise for Raytraced Soft Shadows."_](https://blog.demofox.org/2020/05/16/using-blue-noise-for-raytraced-soft-shadows/)

[Alan Wolfe, _"Not All Blue Noise is Created Equal."_](https://blog.demofox.org/2018/08/12/not-all-blue-noise-is-created-equal/)


[^1]: R. Ulichney, _"The void-and-cluster method
for dither array generation"_ (1993).

[^2]: W. Purgathofer, R.F. Tobler & M. Geiler, _"Forced random dithering: improved threshold matrices for ordered dithering"_ (1994).
