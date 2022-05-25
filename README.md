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
WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWNNXXXKKOO0kkxxddoollcc;;;,,,'''...............''''',,,,;;;;ccccclllloooooddddddxxxxxxxxkkkkkkkkkkkkk00000000000000
NWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWNNNXXKKOO00kkxxddolllcc;;,,,'''.................''''',,,,;;;;;cccclllllooooooddddddxxxxxxxkkkkkkkk0000000000OOOOOOO
XNNNNNNNWWWWWWWWWWWWWWWWWWWWWWWWWWNNNXXKKKOO00kxxddoollcc;;;,,''''......     .......''''',,,,;;;;;ccccclllllooooooddddddxxxxxxkkkkkk0000000OOOOOOOKKKK
KXXXXXXXNNNNNNNWWWWWWWWWWWWWWWWWNNNNXXXKKOO00kkxxddoollcc;;,,,'''.....         ......'''',,,,,;;;;;ccccccllllllooooodddddxxxxxkkkkk000000OOOOOKKKKKKXX
OOKKKKKKXXXXXXNNNNNNNWWWWWWWWNNNNNXXXXKKOO000kkxxdoollccc;;,,''''.....         ......''''',,,,,;;;;;;cccccclllllloooooddddxxxxxkkkk0000OOOOOKKKKXXXXXN
00OOOOOOKKKKKXXXXXXNNNNNNNNNNNNNXXXXKKKOOO00kkxxddoollcc;;,,,'''.....          ......''''',,,,,,;;;;;;cccccccllllllooooddddxxxxkkkk000OOOOKKKKXXXXNNNN
k000000OOOOOOKKKKKXXXXXXXXXXXXXXXXXKKKOOO00kkkxxddoollcc;;,,,'''.....          ......''''',,,,,,;;;;;;;;ccccccclllllooooddddxxxkkkk000OOOKKKXXXXNNNWWW
xkkkkkk00000OOOOOKKKKKXXXXXXXXXXXKKKKOOO000kkxxddoollcc;;;,,''''....           ......''''',,,,,,,;;;;;;;;;cccccccllllloooddddxxxkkk000OOKKKXXXNNNWWWWW
xxxxxxkkkkk00000OOOOKKKKKKKKKKKKKKKOOOO000kkxxdddoollcc;;,,,'''.....          ......'''''',,,,,,,,;;;;;;;;;;;ccccccllllooodddxxxkkk00OOOKKXXXNNWWWWWWW
ddddddxxxxkkkkk0000OOOOOOKKKKKKKOOOOOO000kkxxxddoollccc;;,,,'''.....         ......'''''',,,,,,,,,,;;;;;;;;;;;;cccccllllooodddxxkkk00OOKKKXXNNWWWWWWWW

WWWWWWWWWWWWWWWWWWWNXXXXXXXXXKOO000k000000kddoddddxkk00000OOOOO000kxdddddxk000kkxxxxdddxxkk00OKKXXXXXXXXKKXXXNWWWWWWNNXXKOOOOOOOOOOOOO00kkk0000OKKXNWW
XXNNNNNNNNNNNNWWWWWNXXKKKKKKKKO0000000000kkddddddddxk000OOKKKOO0000kxxxxxk00000kxxxxxddxxxk00OOKKKOOOOOOOKKKXXNWWWWWNXKKOOO0000OOKOOO000kkkk00OKXNNWWW
OOKXXXXXXNNNXXNNNNNNXKKKKOOOOO0000000O000kkxdddddddxkk00OKKKKOO0O00kkkkkk00O000kkxxxxxddxxxk00000000kkk00OOOKXNWWNNXXXKKOO000000OOOOO000kkk000OKNNWWWW
000OOOOKKXXXXXXXXXXXKKKKO00000000000OO0000kxdddddddxxkk0OOKKKOOOO00000kk0OOO00kkkkkxxxdddddxkkkkkkkkxxxxk0OOOXNNNXXXXXKKOO0000k0OOOOOO000000OOOXNWWWWW
kkkkkk0OKKKKKKKKKKKKKOOOO0kkkkkk0000OOO00kkxxxxddddxxkk0OOOKKKKOO0000000OKKO00k0000kxxxddoodxxxxxxxxxxddxk00OKXNNXXKKXXKOO000kk00OOOOOOOOOOOOKKXNWWWWW
dddddxk0OOOOOOO0000OOOOO00kkkxxk00000OOO0kkxxxxddddxxkk0OOKKKKKKOO00000OOKKKO0000000kxxdooooddxxxxxxdddddxk0OKXXNXXKKKKKOO000kkk000OKKKOKKKKKKXNWWWWWW
ddooodxk0OO0OO0kkkk000000k00kxxkk000OKKO0kkxkkkxddddxxk00OKXXXKKKO0000OOKKKKKO000000kkxdollloodddddddoooodxk0KKXXXKKOOOOOO0000kkk0OKXXKKKKKKKKXNWWWWWW
oooloodx000000kkxxkk00000000kxxkkk0OKKXKO0kkk0kkxxxddxkk0OKXXXXKKKO000OOKXKKOO0000000kxdollllloodddooooooodx0OKKKKOOOOOOOO00000000OKXXXXXKKKKXXNWWWWWW
oolllooxk000kkkxxxxkkkk00000kkxxkk0KXXXKO000000kkkxxxkk00OXXNXXXXKKO000OKXKKOO0000000kkdoooollllooooolllooodkOOOOOOOOOOOOO00OOO00OOKKKKXXKKKKXNWWWWWWW
ooooooodxkkkkkxxxxxkkkk00000kkkkkkOKXNXXKOOO00000kkxk0000OKXNNXXXKKO000OKXKOOOOOOO0000kxdddolllllooollllllodx00OOO000OOO0000OOOOOOOKKOOKKKKKXXNWWWWWWW

WWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWNNNXXKKOOO00kkxxddoollcc;;;,,,,'''''........'''''',,,,;;;;ccccclllloooooddddddxxxxxxxxkkkkkkkkkkkkk00000000000000
NWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWWNNNXXKKKOO00kkxxddoollccc;;;,,,'''''........''''''',,,,;;;;;cccclllllooooooddddddxxxxxxxkkkkkkkk0000000000OOOOOOO
XNNNNNNNWWWWWWWWWWWWWWWWWWWWWWWWWWWWNNNNXXKKKOO00kkxxxddoollccc;;,,,,'''''........'''''',,,,,;;;;;ccccclllllooooooddddddxxxxxxkkkkkk0000000OOOOOOOKKKK
KXXXXXXXNNNNNNNWWWWWWWWWWWWWWWWWWWWWNNNNXXXKKOOO00kkxxddooollcc;;;,,,,''''''.....''''''',,,,,,;;;;;ccccccllllllooooodddddxxxxxkkkkk000000OOOOOKKKKKKXX
OOKKKKKKXXXXXXNNNNNNNWWWWWWWWWWWWWWNNNNNXXXKKKOOO00kkxxddooollccc;;,,,,''''''''''''''''',,,,,,,;;;;;;cccccclllllloooooddddxxxxxkkkk0000OOOOOKKKKXXXXXN
00OOOOOOKKKKKXXXXXXNNNNNNNWWWWWWWWWNNNNNNXXXKKKOOO00kkxxddooollccc;;;,,,,''''''''''''''',,,,,,,;;;;;;;cccccccllllllooooddddxxxxkkkk000OOOOKKKKXXXXNNNN
k000000OOOOOOKKKKKXXXXXNNNNNNNNNNNNNNNNNNNXXXKKKOOO00kkxxxddoolllcc;;;;,,,,,''''''''''',,,,,,,,,;;;;;;;;ccccccclllllooooddddxxxkkkk000OOOKKKXXXXNNNWWW
xkkkkkk00000OOOOKKKKKXXXXXNNNNNNNNNNNNNNNNNXXXXKKKOO00kkkxxddooollccc;;;;,,,,,,,,,,,,,,,,,,,,,,,;;;;;;;;;;cccccccllllloooddddxxxkkk000OOKKKXXXNNNWWWWW
xxxxxxkkkkk0000OOOOKKKKXXXXXNNNNNNNNNNNNNNNNNXXXKKKOOO00kkxxxddoolllccc;;;;;,,,,,,,,,,,,,,,,,,,,,;;;;;;;;;;;;ccccccllllooodddxxxkkk00OOOKKXXXNNWWWWWWW
ddddddxxxxkkkk0000OOOKKKKXXXXNNNNNNWWWWWWWNNNNNXXXKKKOOO00kkxxdddoolllcccc;;;;;;;,,,,,,,,,,,,,,,,;;;;;;;;;;;;;;cccccllllooodddxxkkk00OOKKKXXNNWWWWWWWW
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
pub fn new_noise(c FastNoiseConfig) // returns noise struct with settings to generate noise

pub fn (fast FastNoiseLite) get_noise_2(xx f64, yy f64) f64 // returns 2D noise in range [-1, 1]

pub fn (fast FastNoiseLite) get_noise_3(xx f64, yy f64, zz f64) f64 // returns 3D noise in range [-1, 1]

pub fn (fast FastNoiseLite) domain_warp_2(mut x &f64, mut y &f64)

pub fn (fast FastNoiseLite) domain_warp_3(mut x &f64, mut y &f64, mut z &f64)

pub fn text_noise(c FastNoiseConfig, width int, height int, warp bool) // see the noise as text
```
