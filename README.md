# fast_noise_lite
A V implementation of  https://github.com/Auburn/FastNoiseLite

Do you want noise? Do you want it fast? Do you want it lite? Well you found it.

## Install
With VPM

`v install depperm.fast_noise_lite`

## Usage

First import the module

`import depperm.fast_noise_lite as fnl`

Then get noise generator

`mut fast := fnl.new_noise() // can pass config parameters`

Finally get the type of noise you'd like (2D or 3D ... no 4D, why would you want that)

`n := fast.get_noise_2(x, y) // returns a float in range [-1, 1]`

## Examples

If you want quick visual representation of various settings use the `text_noise` function, with config, dimensions, and warping or not:

`fnl.text_noise(fnl.FastNoiseConfig{}, 150, 50, false)`

TODO

## Docs

TODO
