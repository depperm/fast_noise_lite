# fast_noise_lite
A V implementation of  https://github.com/Auburn/FastNoiseLite

Do you want noise? Do you want it fast? Do you want it lite? Well you found it.

## Install
With VPM

`v install depperm.noise`

## Usage

First import the module

`import depperm.noise as fnl`

Then get noise generator

`mut fast := fnl.new_noise() // can pass config parameters`

Finally get the type of noise you'd like (2D or 3D ... no 4D, why would you want that)

`n := fast.get_noise_2(x, y) // returns a float in range [-1, 1]`

## Examples

If you want quick visual representation of various settings use the `text_noise` function, with config (see docs), dimensions (width, height), and warping or not:

`fnl.text_noise(fnl.FastNoiseConfig{}, 150, 50, false)`

OR

```
import depperm.noise as fnl

fn main(){
    
    fnl.text_noise(fnl.FastNoiseConfig{},150,10,false)
    
    println('')
    fnl.text_noise(fnl.FastNoiseConfig{
        m_fractal_type: fnl.FractalType.fbm
        m_octaves: 5
    },150,10,false)

    println('')
    fnl.text_noise(fnl.FastNoiseConfig{
        m_seed: 101
    },150,10,false)
}
```
Which outputs:
```
ddddddddxxxxxxxxxxxxxxxxxxxxxxddddddooooollllccc;;;;,,,,'''''......               .........''''''''',,,,,,,,,;;;;;;;;;;;;;;ccccccccccccccccccccccccccc
oddddddddddddddddxxxxxxxxxdddddddddooooollllcccc;;;;,,,,''''......                 .........''''''''',,,,,,,,,,,;;;;;;;;;;;;;cccccccccccccccccclllllll
ooooooooddddddddddddddddddddddddddooooolllllccc;;;;,,,,'''''......                  .........'''''''''',,,,,,,,,,,;;;;;;;;;;;;ccccccccccccclllllllllll
loooooooooooooodddddddddddddddddooooooollllcccc;;;;,,,,''''......                    .........''''''''''',,,,,,,,,,,;;;;;;;;;;cccccccccccllllllllllloo
lllllllloooooooooooooddddddddooooooooollllccccc;;;,,,,'''''......                    ..........'''''''''''',,,,,,,,,,,;;;;;;;;;ccccccccllllllllloooooo
cclllllllllllooooooooooooooooooooooollllllcccc;;;;,,,,''''......                     ...........''''''''''''',,,,,,,,,,;;;;;;;;ccccccclllllllloooooooo
ccccccclllllllllllooooooooooooooooollllllccccc;;;;,,,,''''......                     ...........''''''''''''''',,,,,,,,,;;;;;;;ccccccclllllloooooooddd
;cccccccccccllllllllllooooooooooolllllllccccc;;;;,,,,'''''......                     ............'''''''''''''''',,,,,,,,;;;;;;;cccccclllllooooooddddd
;;;;;;cccccccccclllllllllllllllllllllllccccc;;;;;,,,,''''......                     ..............''''''''''''''''',,,,,,,;;;;;;cccccllllloooooddddddd
;;;;;;;;;;ccccccccclllllllllllllllllllccccc;;;;;,,,,'''''......                    ................''''''''''''''''',,,,,,,;;;;;cccccllllloooodddddxxx

dddddddddddddddddddooooooooooollccccccccccc;;,;;;;;ccccccclllllcccc;;;;;;;cccccc;;;;;;;;;ccccllloooooooolooooodddxdddooollllllllllllllcccccccccllooodd
oooooooooodooodddddooolllllllllcccccccccccc;;;;;;;;;cccclllllllccccc;;;;;ccccccc;;;;;;;;;;ccclllllllllllllllooodddddoolllllclcclllllllcccccccclloooddd
llloooooooooooooooooolllllllllccccccclccccc;;;;;;;;;cccclllllllcllccccccccllccccc;;;;;;;;;;ccccccccccccccllllooddooooollllcccccclllllccccccccclloodddx
cccllllllooooooooooolllllcccccccccccllccccc;;;;;;;;;;ccclllllllllcccccccclllcccccccc;;;;;;;;ccccccccc;;;ccllloooooooooolllccccccllllllccccccllloodddxx
;cccccclllollllllllllllllccccccccccclllcccc;;;;;;;;;;ccclllllllllccccccclllllccccccc;;;;;,,;;;;;;;;;;;;;;cccllooooolloolllcccccccllllllllllllllooddxxx
;;;;;;cclllllllcccllllllccccc;cccccclllcccc;;;;;;;;;;cccclllllllllcccclllllllcccccccc;;;,,,,;;;;;;;;;;;;;;ccllooooolllllllccccccccclllllllllllooddxxxx
;;,,,;;cclllllccccccccccccccc;;cccccllllccc;ccc;;;;;;;ccclloooolllccccllllllllcccccccc;;,,,,,,;;;;;;;;,,;;;ccllooolllllllllcccccccclooolllllllooddxxxk
,,,,,,;;cccccccc;;ccccccccccc;;cccclllllcccccccc;;;;;;ccclloooollllccclllollllcccccccc;;,,,,,,,,;;;,,,,,,,;;cllllllllllllllcccccccllooooolllloooddxxkk
,,,,,,,;ccccccc;;;;ccccccccccc;;ccllooollccccccccc;;;cccclloooooolllcccllollllccccccccc;,,,,,,,,,,,,,,,,,,,;clllllllllllllcclllcclllllllollllooddxxkkk
,,,,,,,;;ccccc;;;;;cccccccccccccccllooolllllccccccc;ccccclloooooolllcccllollllllllccccc;;;;,,,,,,,,,,,,,,,,;;cclllccclllcccclllllllllllllllloooddxxkkk

ddddddddxxxxxxxxxxxxxxxxxxxxxxdddddddooooolllllcccc;;;;,,,,'''''.........        ..........''''''''',,,,,,,,,;;;;;;;;;;;;;;ccccccccccccccccccccccccccc
oddddddddddddddddxxxxxxxxxxddddddddddooooolllllcccc;;;;,,,,''''''........        ...........''''''''',,,,,,,,,,,;;;;;;;;;;;;;cccccccccccccccccclllllll
ooooooooddddddddddddddddddddddddddddoooooolllllcccc;;;;;,,,,'''''.........        ...........'''''''''',,,,,,,,,,,;;;;;;;;;;;;ccccccccccccclllllllllll
loooooooooooooodddddddddddddddddddddooooooolllllcccc;;;;,,,,,'''''..........     .............''''''''''',,,,,,,,,,,;;;;;;;;;;cccccccccccllllllllllloo
lllllllloooooooooooooddddddddddddddoooooooollllllcccc;;;;,,,,,'''''............................'''''''''''',,,,,,,,,,,;;;;;;;;;ccccccccllllllllloooooo
cclllllllllllooooooooooooodddddddddooooooooollllllcccc;;;;,,,,,''''''..........................'''''''''''''',,,,,,,,,,;;;;;;;;ccccccclllllllloooooooo
ccccccclllllllllllooooooooooooooooooooooooooollllllcccc;;;;;,,,,,''''''.........................''''''''''''''',,,,,,,,,;;;;;;;ccccccclllllloooooooddd
;cccccccccccllllllllloooooooooooooooooooooooooolllllccccc;;;;,,,,,'''''''.......................''''''''''''''''',,,,,,,,;;;;;;;cccccclllllooooooddddd
;;;;;;cccccccccllllllllooooooooooooooooooooooooollllllcccc;;;;;,,,,,''''''''.....................'''''''''''''''''',,,,,,,;;;;;;cccccllllloooooddddddd
;;;;;;;;;;ccccccccllllllloooooooooodddddddoooooooollllllcccc;;;;;,,,,,'''''''''''................''''''''''''''''''',,,,,,,;;;;;cccccllllloooodddddxxx
```

## Docs

Inital Config

```
pub struct FastNoiseConfig {
	m_seed             int             = 1337
	m_freqency         f64             = 0.01
	m_noise_type       NoiseType       = NoiseType.open_simplex2
	m_rotation_type3d  RotationType3D  = RotationType3D.@none
	m_transform_type3d TransformType3D = TransformType3D.default_open_simplex2

	m_fractal_type       FractalType = FractalType.@none
	m_octaves            int = 3
	m_lacunarity         f64 = 2.0
	m_gain               f64 = 0.5
	m_weighted_strength  f64 = 0.0
	m_ping_pong_strength f64 = 2.0

	m_fractal_bounding f64 = 1 / 1.75

	m_cellular_distance_function CellularDistanceFunction = CellularDistanceFunction.euclidean_sq
	m_cellular_return_type       CellularReturnType       = CellularReturnType.distance
	m_cullular_jitter_modifier   f64 = 1.0

	m_domain_warp_type      DomainWarpType  = DomainWarpType.open_simplex2
	m_warp_transform_type3d TransformType3D = TransformType3D.default_open_simplex2
	m_domain_warp_amp       f64 = 1.0
	
	m_random_num_range bool
}
```

### Conig enums

```
pub enum NoiseType {
	open_simplex2
	open_simplex2s
	cellular
	perlin
	value_cubic
	value
}

pub enum RotationType3D {
	@none
	improve_xy_planes
	improve_zx_planes
}

pub enum TransformType3D {
	@none
	improve_xy_planes
	improve_zx_planes
	default_open_simplex2
}

pub enum FractalType {
	@none
	fbm
	ridged
	ping_pong
	domain_warp_progressive
	domain_warp_independent
}

pub enum CellularDistanceFunction {
	euclidean
	euclidean_sq
	manhattan
	hybrid
}

pub enum CellularReturnType {
	cell_value
	distance
	distance2
	distance2_add
	distance2_sub
	distance2_mul
	distance2_div
}

pub enum DomainWarpType {
	open_simplex2
	open_simplex2_reduced
	basic_grid
}
```

Public functions

```
pub fn new_noise(c FastNoiseConfig) FastNoiseLite // returns noise struct with settings to generate noise

pub fn (fast FastNoiseLite) get_noise_2(xx f64, yy f64) f64 // returns 2D noise in range [-1, 1]

pub fn (fast FastNoiseLite) get_noise_3(xx f64, yy f64, zz f64) f64 // returns 3D noise in range [-1, 1]

pub fn (fast FastNoiseLite) domain_warp_2(mut x f64, mut y f64)

pub fn (fast FastNoiseLite) domain_warp_3(mut x f64, mut y f64, mut z f64)

pub fn text_noise(c FastNoiseConfig, width int, height int, warp bool) // see the noise as text

// CONFIG individual setters
pub fn (mut fast FastNoiseLite) set_seed(s int)

pub fn (mut fast FastNoiseLite) set_frequency(f f64)

pub fn (mut fast FastNoiseLite) set_noise_type(nt NoiseType)

pub fn (mut fast FastNoiseLite) set_rotation_type(rt RotationType3D)

pub fn (mut fast FastNoiseLite) set_transform_type(tt TransformType3D)

pub fn (mut fast FastNoiseLite) set_fractal_type(ft FractalType)

pub fn (mut fast FastNoiseLite) set_octaves(o int)

pub fn (mut fast FastNoiseLite) set_lacunarity(l f64) 

pub fn (mut fast FastNoiseLite) set_gain(g f64)

pub fn (mut fast FastNoiseLite) set_weighted_strength(ws f64)

pub fn (mut fast FastNoiseLite) set_ping_pong_strength(pps f64)

pub fn (mut fast FastNoiseLite) set_fractal_bounding(fb f64)

pub fn (mut fast FastNoiseLite) set_cellular_distance_function(cdf CellularDistanceFunction)

pub fn (mut fast FastNoiseLite) set_cellular_return_type(crt CellularReturnType)

pub fn (mut fast FastNoiseLite) set_cellular_jitter_modifier(cjm f64)

pub fn (mut fast FastNoiseLite) set_domain_warp_type(dwt DomainWarpType)

pub fn (mut fast FastNoiseLite) set_warp_transform_type(wtt TransformType3D)

pub fn (mut fast FastNoiseLite) set_domain_warp_amp(dwa f64)
```
