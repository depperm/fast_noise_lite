module noise

// MIT License
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files(the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and / or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions :
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.
//
// VERSION: 0.0.1
// source: https://github.com/Auburn/FastNoise
import math

const (
	prime_x = 501125321
	prime_y = 1136930381
	prime_z = 1720413743

	sqrt3   = 1.7320508075688772935274463415059
	g2      = (3.0 - sqrt3) / 6
	f2      = .5 * (sqrt3 - 1)
	r3      = 2.0 / 3
)

fn hash_2(seed int, x_primed int, y_primed int) int {
	mut hash := math.powi(math.powi(seed, x_primed), y_primed)

	hash *= 0x27d4eb2d
	return int(hash)
}

fn hash_3(seed int, x_primed int, y_primed int, z_primed int) int { // TODO change to i64?
	mut hash := math.powi(math.powi(math.powi(seed, x_primed), y_primed), z_primed)

	hash *= 0x27d4eb2d
	return int(hash)
}

fn val_coord_2(seed int, x_primed int, y_primed int) f64 {
	mut hash := u32(hash_2(seed, x_primed, y_primed))

	hash *= hash
	hash ^= hash << 19
	return f64(hash) * (1.0 / 2147483648.0)
}

fn val_coord_3(seed int, x_primed int, y_primed int, z_primed int) f64 {
	mut hash := u32(hash_3(seed, x_primed, y_primed, z_primed))

	hash *= hash
	hash ^= hash << 19
	return f64(hash) * (1.0 / 2147483648.0)
}

fn grad_coord_2(seed int, x_primed int, y_primed int, xd f64, yd f64) f64 {
	mut hash := u32(hash_2(seed, x_primed, y_primed))
	hash ^= hash >> 15
	hash &= 127 << 1

	xg := gradients2d[hash]
	yg := gradients2d[hash | 1]

	return xd * xg + yd * yg
}

fn grad_coord_3(seed int, x_primed int, y_primed int, z_primed int, xd f64, yd f64, zd f64) f64 {
	mut hash := u32(hash_3(seed, x_primed, y_primed, z_primed))
	hash ^= hash >> 15
	hash &= 63 << 2

	xg := gradients3d[hash]
	yg := gradients3d[hash | 1]
	zg := gradients3d[hash | 2]

	return xd * xg + yd * yg + zd * zg
}

fn grad_coord_out_2(seed int, x_primed int, y_primed int, mut xo &f64, mut yo &f64) {
	hash := u32(hash_2(seed, x_primed, y_primed)) & (255 << 1)

	xo = rand_vecs2d[hash]
	yo = rand_vecs2d[hash | 1]
}

fn grad_coord_out_3(seed int, x_primed int, y_primed int, z_primed int, mut xo &f64, mut yo &f64, mut zo &f64) {
	hash := u32(hash_3(seed, x_primed, y_primed, z_primed)) & (255 << 2)

	xo = rand_vecs3d[hash]
	yo = rand_vecs3d[hash | 1]
	zo = rand_vecs3d[hash | 2]
}

fn grad_coord_dual_2(seed int, x_primed int, y_primed int, xd f64, yd f64, mut xo &f64, mut yo &f64) {
	hash := u32(hash_2(seed, x_primed, y_primed))
	index1 := hash & (127 << 1)
	index2 := (hash >> 7) & (255 << 1)

	xg := gradients2d[index1]
	yg := gradients2d[index1 | 1]
	value := xd * xg + yd * yg

	xgo := rand_vecs2d[index2]
	ygo := rand_vecs2d[index2 | 1]

	xo = value * xgo
	yo = value * ygo
}

fn grad_coord_dual_3(seed int, x_primed int, y_primed int, z_primed int, xd f64, yd f64, zd f64, mut xo &f64, mut yo &f64, mut zo &f64) {
	hash := u32(hash_3(seed, x_primed, y_primed, z_primed))
	index1 := hash & (63 << 2)
	index2 := (hash >> 6) & (255 << 2)

	xg := gradients3d[index1]
	yg := gradients3d[index1 | 1]
	zg := gradients3d[index1 | 2]
	value := xd * xg + yd * yg + zd * zg

	xgo := rand_vecs3d[index2]
	ygo := rand_vecs3d[index2 | 1]
	zgo := rand_vecs3d[index2 | 2]

	xo = value * xgo
	yo = value * ygo
	zo = value * zgo
}

fn ping_pong(t f64) f64 {
	tt := t - int(t * .5) * 2
	return ternary<f64>(tt < 1, tt, 2 - tt)
}

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

[params]
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
	m_cellular_jitter_modifier   f64 = 1.0

	m_domain_warp_type      DomainWarpType  = DomainWarpType.open_simplex2
	m_warp_transform_type3d TransformType3D = TransformType3D.default_open_simplex2
	m_domain_warp_amp       f64 = 1.0

	m_random_num_range bool = false
}

struct FastNoiseLite {
mut:
	m_seed             int
	m_freqency         f64
	m_noise_type       NoiseType
	m_rotation_type3d  RotationType3D
	m_transform_type3d TransformType3D

	m_fractal_type       FractalType
	m_octaves            int
	m_lacunarity         f64
	m_gain               f64
	m_weighted_strength  f64
	m_ping_pong_strength f64

	m_fractal_bounding f64

	m_cellular_distance_function CellularDistanceFunction
	m_cellular_return_type       CellularReturnType
	m_cellular_jitter_modifier   f64

	m_domain_warp_type      DomainWarpType
	m_warp_transform_type3d TransformType3D
	m_domain_warp_amp       f64

	m_random_num_range bool
}

pub fn (mut fast FastNoiseLite) set_seed(s int) {
	fast.m_seed = s
}

pub fn (mut fast FastNoiseLite) set_frequency(f f64) {
	fast.m_freqency = f
}

pub fn (mut fast FastNoiseLite) set_noise_type(nt NoiseType) {
	fast.m_noise_type = nt
}

pub fn (mut fast FastNoiseLite) set_rotation_type(rt RotationType3D) {
	fast.m_rotation_type3d = rt
}

pub fn (mut fast FastNoiseLite) set_transform_type(tt TransformType3D) {
	fast.m_transform_type3d = tt
}

pub fn (mut fast FastNoiseLite) set_fractal_type(ft FractalType) {
	fast.m_fractal_type = ft
}

pub fn (mut fast FastNoiseLite) set_octaves(o int) {
	fast.m_octaves = o
}

pub fn (mut fast FastNoiseLite) set_lacunarity(l f64) {
	fast.m_lacunarity = l
}

pub fn (mut fast FastNoiseLite) set_gain(g f64) {
	fast.m_gain = g
}

pub fn (mut fast FastNoiseLite) set_weighted_strength(ws f64) {
	fast.m_weighted_strength = ws
}

pub fn (mut fast FastNoiseLite) set_ping_pong_strength(pps f64) {
	fast.m_ping_pong_strength = pps
}

pub fn (mut fast FastNoiseLite) set_fractal_bounding(fb f64) {
	fast.m_fractal_bounding = fb
}

pub fn (mut fast FastNoiseLite) set_cellular_distance_function(cdf CellularDistanceFunction) {
	fast.m_cellular_distance_function = cdf
}

pub fn (mut fast FastNoiseLite) set_cellular_return_type(crt CellularReturnType) {
	fast.m_cellular_return_type = crt
}

pub fn (mut fast FastNoiseLite) set_cellular_jitter_modifier(cjm f64) {
	fast.m_cellular_jitter_modifier = cjm
}

pub fn (mut fast FastNoiseLite) set_domain_warp_type(dwt DomainWarpType) {
	fast.m_domain_warp_type = dwt
}

pub fn (mut fast FastNoiseLite) set_warp_transform_type(wtt TransformType3D) {
	fast.m_warp_transform_type3d = wtt
}

pub fn (mut fast FastNoiseLite) set_domain_warp_amp(dwa f64) {
	fast.m_domain_warp_amp = dwa
}

// new_noise returns FastNoiseLite struct to generate noise
//
// Example:
// ```
// mut fast := new_noise()
// n := fast.get_noise_2(x,y)
// ```
pub fn new_noise(c FastNoiseConfig) &FastNoiseLite {
	return &FastNoiseLite{
		m_seed: c.m_seed
		m_freqency: c.m_freqency
		m_noise_type: c.m_noise_type
		m_rotation_type3d: c.m_rotation_type3d
		m_transform_type3d: c.m_transform_type3d
		m_fractal_type: c.m_fractal_type
		m_octaves: c.m_octaves
		m_lacunarity: c.m_lacunarity
		m_gain: c.m_gain
		m_weighted_strength: c.m_weighted_strength
		m_ping_pong_strength: c.m_ping_pong_strength
		m_fractal_bounding: c.m_fractal_bounding
		m_cellular_distance_function: c.m_cellular_distance_function
		m_cellular_return_type: c.m_cellular_return_type
		m_cellular_jitter_modifier: c.m_cellular_jitter_modifier
		m_domain_warp_type: c.m_domain_warp_type
		m_warp_transform_type3d: c.m_warp_transform_type3d
		m_domain_warp_amp: c.m_domain_warp_amp
	}
}

// Noise

// get_noise_2 calculates 2D noise at given position using current settings
//
// returns [-1,1]
pub fn (fast FastNoiseLite) get_noise_2(xx f64, yy f64) f64 {
	mut x, mut y := xx, yy
	fast.transform_noise_coord_2(mut x, mut y)

	mut n := 0.0
	match fast.m_fractal_type {
		.fbm {
			n = fast.gen_fractal_fbm_2(x, y)
		}
		.ridged {
			n = fast.gen_fractal_ridged_2(x, y)
		}
		.ping_pong {
			n = fast.gen_fractal_ping_pong_2(x, y)
		}
		else {
			n = fast.gen_noise_single_2(fast.m_seed, x, y)
		}
	}
	if fast.m_random_num_range {
		n += 1
		n /= 2
	}
	return n
}

// get_noise_3 calculates 3D noise at given position using current settings
//
// returns [-1,1]
pub fn (fast FastNoiseLite) get_noise_3(xx f64, yy f64, zz f64) f64 {
	mut x, mut y, mut z := xx, yy, zz
	fast.transform_noise_coord_3(mut x, mut y, mut z)

	mut n := 0.0
	match fast.m_fractal_type {
		.fbm {
			n = fast.gen_fractal_fbm_3(x, y, z)
		}
		.ridged {
			n = fast.gen_fractal_ridged_3(x, y, z)
		}
		.ping_pong {
			n = fast.gen_fractal_ping_pong_3(x, y, z)
		}
		else {
			n = fast.gen_noise_single_3(fast.m_seed, x, y, z)
		}
	}
	if fast.m_random_num_range {
		n += 1
		n /= 2
	}
	return n
}

fn (fast FastNoiseLite) transform_noise_coord_2(mut x &f64, mut y &f64) {
	x *= fast.m_freqency
	y *= fast.m_freqency

	match fast.m_noise_type {
		.open_simplex2, .open_simplex2s {
			t := (*x + *y) * noise.f2
			x += t
			y += t
		}
		else {}
	}
}

fn (fast FastNoiseLite) transform_noise_coord_3(mut x &f64, mut y &f64, mut z &f64) {
	x *= fast.m_freqency
	y *= fast.m_freqency
	z *= fast.m_freqency

	match fast.m_transform_type3d {
		.improve_xy_planes {
			xy := *x + *y
			s2 := xy * -0.211324865405187
			z *= 0.577350269189626
			x += s2 - z
			y += s2 - z
			z += xy * 0.577350269189626
		}
		.improve_zx_planes {
			xz := *x + *z
			s2 := xz * -0.211324865405187
			y *= 0.577350269189626
			x += s2 - y
			z += s2 - y
			y += xz * 0.577350269189626
		}
		.default_open_simplex2 {
			r := (*x + *y + *z) * noise.r3
			x = r - x
			y = r - y
			z = r - z
		}
		else {}
	}
}

fn (fast FastNoiseLite) gen_noise_single_2(seed int, x f64, y f64) f64 {
	match fast.m_noise_type {
		.open_simplex2 {
			return fast.single_simplex(seed, x, y)
		}
		.open_simplex2s {
			return fast.single_open_simplex2s_2(seed, x, y)
		}
		.cellular {
			return fast.single_cellular_2(seed, x, y)
		}
		.perlin {
			return fast.single_perlin_2(seed, x, y)
		}
		.value_cubic {
			return fast.single_value_cubic_2(seed, x, y)
		}
		.value {
			return fast.single_value_2(seed, x, y)
		}
	}
}

fn (fast FastNoiseLite) gen_noise_single_3(seed int, x f64, y f64, z f64) f64 {
	match fast.m_noise_type {
		.open_simplex2 {
			return fast.single_open_simplex2(seed, x, y, z)
		}
		.open_simplex2s {
			return fast.single_open_simplex2s_3(seed, x, y, z)
		}
		.cellular {
			return fast.single_cellular_3(seed, x, y, z)
		}
		.perlin {
			return fast.single_perlin_3(seed, x, y, z)
		}
		.value_cubic {
			return fast.single_value_cubic_3(seed, x, y, z)
		}
		.value {
			return fast.single_value_3(seed, x, y, z)
		}
	}
}

fn (fast FastNoiseLite) gen_fractal_fbm_2(xx f64, yy f64) f64 {
	mut x, mut y := xx, yy
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := fast.gen_noise_single_2(seed, x, y)
		seed++
		sum += noise * amp
		amp *= lerp(1.0, fast_min(noise + 1, 2) * .5, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		amp *= fast.m_gain
	}
	return sum
}

fn (fast FastNoiseLite) gen_fractal_fbm_3(xx f64, yy f64, zz f64) f64 {
	mut x, mut y, mut z := xx, yy, zz
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := fast.gen_noise_single_3(seed, x, y, z)
		seed++
		sum += noise * amp
		amp *= lerp(1.0, (noise + 1) * .5, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		z *= fast.m_lacunarity
		amp *= fast.m_gain
	}
	return sum
}

fn (fast FastNoiseLite) gen_fractal_ridged_2(xx f64, yy f64) f64 {
	mut x, mut y := xx, yy
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := fast_abs(fast.gen_noise_single_2(seed, x, y))
		seed++
		sum += (noise * -2 + 1) * amp
		amp *= lerp(1.0, 1 - noise, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		amp *= fast.m_gain
	}

	return sum
}

fn (fast FastNoiseLite) gen_fractal_ridged_3(xx f64, yy f64, zz f64) f64 {
	mut x, mut y, mut z := xx, yy, zz
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := fast_abs(fast.gen_noise_single_3(seed, x, y, z))
		seed++
		sum += (noise * -2 + 1) * amp
		amp *= lerp(1.0, 1 - noise, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		z *= fast.m_lacunarity
		amp *= fast.m_gain
	}

	return sum
}

fn (fast FastNoiseLite) gen_fractal_ping_pong_2(xx f64, yy f64) f64 {
	mut x, mut y := xx, yy
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := ping_pong((fast.gen_noise_single_2(seed, x, y) + 1) * fast.m_ping_pong_strength)
		seed++
		sum += (noise - 0.5) * 2 * amp
		amp *= lerp(1.0, noise, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		amp *= fast.m_gain
	}

	return sum
}

fn (fast FastNoiseLite) gen_fractal_ping_pong_3(xx f64, yy f64, zz f64) f64 {
	mut x, mut y, mut z := xx, yy, zz
	mut seed := fast.m_seed
	mut sum := 0.0
	mut amp := fast.m_fractal_bounding

	for _ in 0 .. fast.m_octaves {
		noise := ping_pong((fast.gen_noise_single_3(seed, x, y, z) + 1) * fast.m_ping_pong_strength)
		seed++
		sum += (noise - 0.5) * 2 * amp
		amp *= lerp(1.0, noise, fast.m_weighted_strength)

		x *= fast.m_lacunarity
		y *= fast.m_lacunarity
		z *= fast.m_lacunarity
		amp *= fast.m_gain
	}

	return sum
}

fn (fast FastNoiseLite) single_simplex(seed int, x f64, y f64) f64 {
	mut i := fast_floor(x)
	mut j := fast_floor(y)
	xi := f64(x - i)
	yi := f64(y - j)

	t := (xi + yi) * noise.g2
	x0 := xi - t
	y0 := yi - t

	i *= noise.prime_x
	j *= noise.prime_y

	mut n0, mut n1, mut n2 := f64(0), f64(0), f64(0)

	a := .5 - x0 * x0 - y0 * y0
	if a > 0 {
		n0 = (a * a) * (a * a) * grad_coord_2(seed, i, j, x0, y0)
	}

	c := f64(2 * (1 - 2 * noise.g2) * (1 / noise.g2 - 2)) * t +
		(f64(-2 * (1 - 2 * noise.g2) * (1 - 2 * noise.g2)) + a)
	if c > 0 {
		x2 := x0 + (2 * noise.g2 - 1)
		y2 := y0 + (2 * noise.g2 - 1)
		n2 = (c * c) * (c * c) * grad_coord_2(seed, i + noise.prime_x, j + noise.prime_y,
			x2, y2)
	}

	if y0 > x0 {
		x1 := x0 + noise.g2
		y1 := y0 + (noise.g2 - 1)
		b := .5 - x1 * x1 - y1 * y1
		if b > 0 {
			n1 = (b * b) * (b * b) * grad_coord_2(seed, i, j + noise.prime_y, x1, y1)
		}
	} else {
		x1 := x0 + (noise.g2 - 1)
		y1 := y0 + noise.g2
		b := .5 - x1 * x1 - y1 * y1
		if b > 0 {
			n1 = (b * b) * (b * b) * grad_coord_2(seed, i + noise.prime_x, j, x1, y1)
		}
	}

	return (n0 + n1 + n2) * 99.83685446303647
}

fn (fast FastNoiseLite) single_open_simplex2s_2(seed int, x f64, y f64) f64 {
	// 2D OpenSimplex2S case is a modified 2D simplex noise.
	/*
	* --- Skew moved to TransformNoiseCoordinate method ---
        * const FNfloat F2 = 0.5 * (SQRT3 - 1)
        * FNfloat s = (x + y) * F2
        * x += s y += s
	*/

	mut i := fast_floor(x)
	mut j := fast_floor(y)
	xi := x - i
	yi := y - j

	i *= noise.prime_x
	j *= noise.prime_y
	i1 := i + noise.prime_x
	j1 := j + noise.prime_y

	t := (xi + yi) * noise.g2
	x0 := xi - t
	y0 := yi - t

	a0 := (2.0 / 3.0) - x0 * x0 - y0 * y0
	mut value := (a0 * a0) * (a0 * a0) * grad_coord_2(seed, i, j, x0, y0)

	a1 := f64(2 * (1 - 2 * noise.g2) * (1 / noise.g2 - 2)) * t +
		(f64(-2 * (1 - 2 * noise.g2) * (1 - 2 * noise.g2)) + a0)
	x1 := x0 - f64(1 - 2 * noise.g2)
	y1 := y0 - f64(1 - 2 * noise.g2)
	value += (a1 * a1) * (a1 * a1) * grad_coord_2(seed, i1, j1, x1, y1)

	// Nested conditionals were faster than compact bit logic/arithmetic.
	xmyi := xi - yi
	if t > noise.g2 {
		if xi + xmyi > 1 {
			x2 := x0 + f64(3 * noise.g2 - 2)
			y2 := y0 + f64(3 * noise.g2 - 1)
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i + (noise.prime_x << 1),
					j + noise.prime_y, x2, y2)
			}
		} else {
			x2 := x0 + noise.g2
			y2 := y0 + f64(noise.g2 - 1)
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i, j + noise.prime_y,
					x2, y2)
			}
		}

		if yi - xmyi > 1 {
			x3 := x0 + f64(3 * noise.g2 - 1)
			y3 := y0 + f64(3 * noise.g2 - 2)
			a3 := (2.0 / 3.0) - x3 * x3 - y3 * y3
			if a3 > 0 {
				value += (a3 * a3) * (a3 * a3) * grad_coord_2(seed, i + noise.prime_x,
					j + (noise.prime_y << 1), x3, y3)
			}
		} else {
			x3 := x0 + f64(noise.g2 - 1)
			y3 := y0 + noise.g2
			a3 := (2.0 / 3.0) - x3 * x3 - y3 * y3
			if a3 > 0 {
				value += (a3 * a3) * (a3 * a3) * grad_coord_2(seed, i + noise.prime_x,
					j, x3, y3)
			}
		}
	} else {
		if xi + xmyi < 0 {
			x2 := x0 + f64(1 - noise.g2)
			y2 := y0 - noise.g2
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i - noise.prime_x,
					j, x2, y2)
			}
		} else {
			x2 := x0 + f64(noise.g2 - 1)
			y2 := y0 + noise.g2
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i + noise.prime_x,
					j, x2, y2)
			}
		}

		if yi < xmyi {
			x2 := x0 - noise.g2
			y2 := y0 - f64(noise.g2 - 1)
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i, j - noise.prime_y,
					x2, y2)
			}
		} else {
			x2 := x0 + noise.g2
			y2 := y0 + f64(noise.g2 - 1)
			a2 := (2.0 / 3.0) - x2 * x2 - y2 * y2
			if a2 > 0 {
				value += (a2 * a2) * (a2 * a2) * grad_coord_2(seed, i, j + noise.prime_y,
					x2, y2)
			}
		}
	}

	return value * 18.24196194486065
}

fn (fast FastNoiseLite) single_open_simplex2s_3(seed int, x f64, y f64, z f64) f64 {
	// 3D OpenSimplex2S case uses two offset rotated cube grids.
	/*
	* --- Rotation moved to TransformNoiseCoordinate method ---
    * r3 := (2.0 / 3.0)
    * r := (x + y + z) * r3 // Rotation, not skew
    * x := r - x, y := r - y, z := r - z
	*/

	mut i := fast_floor(x)
	mut j := fast_floor(y)
	mut k := fast_floor(z)
	xi := x - i
	yi := y - j
	zi := z - k

	i *= noise.prime_x
	j *= noise.prime_y
	k *= noise.prime_z
	seed2 := seed + 1293373

	x_n_mask := int(-0.5 - xi)
	y_n_mask := int(-0.5 - yi)
	z_n_mask := int(-0.5 - zi)

	x0 := xi + x_n_mask
	y0 := yi + y_n_mask
	z0 := zi + z_n_mask
	a0 := 0.75 - x0 * x0 - y0 * y0 - z0 * z0
	mut value := (a0 * a0) * (a0 * a0) * grad_coord_3(seed, i + (x_n_mask & noise.prime_x),
		j + (y_n_mask & noise.prime_y), k + (z_n_mask & noise.prime_z), x0, y0, z0)

	x1 := xi - 0.5
	y1 := yi - 0.5
	z1 := zi - 0.5
	a1 := 0.75 - x1 * x1 - y1 * y1 - z1 * z1
	value += (a1 * a1) * (a1 * a1) * grad_coord_3(seed2, i + noise.prime_x, j + noise.prime_y,
		k + noise.prime_z, x1, y1, z1)

	x_a_flip_mask_0 := ((u32(x_n_mask) | 1) << 1) * x1
	y_a_flip_mask_0 := ((u32(y_n_mask) | 1) << 1) * y1
	z_a_flip_mask_0 := ((u32(z_n_mask) | 1) << 1) * z1
	x_a_flip_mask_1 := (-2 - (u32(x_n_mask) << 2)) * x1 - 1.0
	y_a_flip_mask_1 := (-2 - (u32(y_n_mask) << 2)) * y1 - 1.0
	z_a_flip_mask_1 := (-2 - (u32(z_n_mask) << 2)) * z1 - 1.0

	mut skip5 := false
	a2 := x_a_flip_mask_0 + a0
	if a2 > 0 {
		x2 := x0 - (x_n_mask | 1)
		y2 := y0
		z2 := z0
		value += (a2 * a2) * (a2 * a2) * grad_coord_3(seed, i + (~x_n_mask & noise.prime_x),
			j + (y_n_mask & noise.prime_y), k + (z_n_mask & noise.prime_z), x2, y2, z2)
	} else {
		a3 := y_a_flip_mask_0 + z_a_flip_mask_0 + a0
		if a3 > 0 {
			x3 := x0
			y3 := y0 - (y_n_mask | 1)
			z3 := z0 - (z_n_mask | 1)
			value += (a3 * a3) * (a3 * a3) * grad_coord_3(seed, i + (x_n_mask & noise.prime_x),
				j + (~y_n_mask & noise.prime_y), k + (~z_n_mask & noise.prime_z), x3,
				y3, z3)
		}

		a4 := x_a_flip_mask_1 + a1
		if a4 > 0 {
			x4 := (x_n_mask | 1) + x1
			y4 := y1
			z4 := z1
			value += (a4 * a4) * (a4 * a4) * grad_coord_3(seed2, i +
				(x_n_mask & (noise.prime_x * 2)), j + noise.prime_y, k + noise.prime_z,
				x4, y4, z4)
			skip5 = true
		}
	}

	mut skip9 := false
	a6 := y_a_flip_mask_0 + a0
	if a6 > 0 {
		x6 := x0
		y6 := y0 - (y_n_mask | 1)
		z6 := z0
		value += (a6 * a6) * (a6 * a6) * grad_coord_3(seed, i + (x_n_mask & noise.prime_x),
			j + (~y_n_mask & noise.prime_y), k + (z_n_mask & noise.prime_z), x6, y6, z6)
	} else {
		a7 := x_a_flip_mask_0 + z_a_flip_mask_0 + a0
		if a7 > 0 {
			x7 := x0 - (x_n_mask | 1)
			y7 := y0
			z7 := z0 - (z_n_mask | 1)
			value += (a7 * a7) * (a7 * a7) * grad_coord_3(seed, i + (~x_n_mask & noise.prime_x),
				j + (y_n_mask & noise.prime_y), k + (~z_n_mask & noise.prime_z), x7, y7,
				z7)
		}

		a8 := y_a_flip_mask_1 + a1
		if a8 > 0 {
			x8 := x1
			y8 := (y_n_mask | 1) + y1
			z8 := z1
			value += (a8 * a8) * (a8 * a8) * grad_coord_3(seed2, i + noise.prime_x, j +
				(y_n_mask & (noise.prime_y << 1)), k + noise.prime_z, x8, y8, z8)
			skip9 = true
		}
	}

	mut skip_bd := false
	a_ba := z_a_flip_mask_0 + a0
	if a_ba > 0 {
		x_ba := x0
		y_ba := y0
		z_ba := z0 - (z_n_mask | 1)
		value += (a_ba * a_ba) * (a_ba * a_ba) * grad_coord_3(seed, i + (x_n_mask & noise.prime_x),
			j + (y_n_mask & noise.prime_y), k + (~z_n_mask & noise.prime_z), x_ba, y_ba,
			z_ba)
	} else {
		a_bb := x_a_flip_mask_0 + y_a_flip_mask_0 + a0
		if a_bb > 0 {
			x_bb := x0 - (x_n_mask | 1)
			y_bb := y0 - (y_n_mask | 1)
			z_bb := z0
			value += (a_bb * a_bb) * (a_bb * a_bb) * grad_coord_3(seed, i +
				(~x_n_mask & noise.prime_x), j + (~y_n_mask & noise.prime_y), k +
				(z_n_mask & noise.prime_z), x_bb, y_bb, z_bb)
		}

		a_bc := z_a_flip_mask_1 + a1
		if a_bc > 0 {
			x_c := x1
			y_c := y1
			z_c := (z_n_mask | 1) + z1
			value += (a_bc * a_bc) * (a_bc * a_bc) * grad_coord_3(seed2, i + noise.prime_x,
				j + noise.prime_y, k + (z_n_mask & (noise.prime_z << 1)), x_c, y_c, z_c)
			skip_bd = true
		}
	}

	if !skip5 {
		a5 := y_a_flip_mask_1 + z_a_flip_mask_1 + a1
		if a5 > 0 {
			x5 := x1
			y5 := (y_n_mask | 1) + y1
			z5 := (z_n_mask | 1) + z1
			value += (a5 * a5) * (a5 * a5) * grad_coord_3(seed2, i + noise.prime_x, j +
				(y_n_mask & (noise.prime_y << 1)), k + (z_n_mask & (noise.prime_z << 1)),
				x5, y5, z5)
		}
	}

	if !skip9 {
		a9 := x_a_flip_mask_1 + z_a_flip_mask_1 + a1
		if a9 > 0 {
			x9 := (x_n_mask | 1) + x1
			y9 := y1
			z9 := (z_n_mask | 1) + z1
			value += (a9 * a9) * (a9 * a9) * grad_coord_3(seed2, i +
				(x_n_mask & (noise.prime_x * 2)), j + noise.prime_y, k +
				(z_n_mask & (noise.prime_z << 1)), x9, y9, z9)
		}
	}

	if !skip_bd {
		a_bd := x_a_flip_mask_1 + y_a_flip_mask_1 + a1
		if a_bd > 0 {
			x_d := (x_n_mask | 1) + x1
			y_d := (y_n_mask | 1) + y1
			z_d := z1
			value += (a_bd * a_bd) * (a_bd * a_bd) * grad_coord_3(seed2, i +
				(x_n_mask & (noise.prime_x << 1)), j + (y_n_mask & (noise.prime_y << 1)),
				k + noise.prime_z, x_d, y_d, z_d)
		}
	}

	return value * 9.046026385208288
}

fn (fast FastNoiseLite) single_open_simplex2(s int, x f64, y f64, z f64) f64 {
	// 3D OpenSimplex2 case uses two offset rotated cube grids.
	/*
	* --- Rotation moved to TransformNoiseCoordinate method ---
        * const FNfloat R3 = (FNfloat)(2.0 / 3.0)
        * FNfloat r = (x + y + z) * R3; // Rotation, not skew
        * x = r - x; y = r - y; z = r - z
	*/
	mut seed := s
	mut i := fast_round(x)
	mut j := fast_round(y)
	mut k := fast_round(z)
	mut x0 := x - i
	mut y0 := y - j
	mut z0 := z - k

	mut x_n_sign := int(-1.0 - x0) | 1
	mut y_n_sign := int(-1.0 - y0) | 1
	mut z_n_sign := int(-1.0 - z0) | 1

	mut ax0 := x_n_sign * -x0
	mut ay0 := y_n_sign * -y0
	mut az0 := z_n_sign * -z0

	i *= noise.prime_x
	j *= noise.prime_y
	k *= noise.prime_z

	mut value := 0.0
	mut a := (0.6 - x0 * x0) - (y0 * y0 + z0 * z0)

	for l := 0; true; l++ {
		if a > 0 {
			value += (a * a) * (a * a) * grad_coord_3(seed, i, j, k, x0, y0, z0)
		}

		mut b := a + 1
		mut i1 := i
		mut j1 := j
		mut k1 := k
		mut x1 := x0
		mut y1 := y0
		mut z1 := z0

		if ax0 >= ay0 && ax0 >= az0 {
			x1 += x_n_sign
			b -= x_n_sign * 2 * x1
			i1 -= x_n_sign * noise.prime_x
		} else if ay0 > ax0 && ay0 >= az0 {
			y1 += y_n_sign
			b -= y_n_sign * 2 * y1
			j1 -= y_n_sign * noise.prime_y
		} else {
			z1 += z_n_sign
			b -= z_n_sign * 2 * z1
			k1 -= z_n_sign * noise.prime_z
		}

		if b > 0 {
			value += (b * b) * (b * b) * grad_coord_3(seed, i1, j1, k1, x1, y1, z1)
		}

		if l == 1 {
			break
		}

		ax0 = 0.5 - ax0
		ay0 = 0.5 - ay0
		az0 = 0.5 - az0

		x0 = x_n_sign * ax0
		y0 = y_n_sign * ay0
		z0 = z_n_sign * az0

		a += (0.75 - ax0) - (ay0 + az0)

		i += (x_n_sign >> 1) & noise.prime_x
		j += (y_n_sign >> 1) & noise.prime_y
		k += (z_n_sign >> 1) & noise.prime_z

		x_n_sign = -x_n_sign
		y_n_sign = -y_n_sign
		z_n_sign = -z_n_sign

		seed = ~seed
	}

	return value * 32.69428253173828125
}

fn (fast FastNoiseLite) single_cellular_2(seed int, x f64, y f64) f64 {
	xr := fast_round(x)
	yr := fast_round(y)

	mut distance0 := 1e10 // f64
	mut distance1 := 1e10
	mut closest_hash := 0

	cellular_jitter := 0.43701595 * fast.m_cellular_jitter_modifier

	mut x_primed := int((xr - 1) * noise.prime_x)
	y_primed_base := int((yr - 1) * noise.prime_y)

	match fast.m_cellular_distance_function {
		.manhattan {
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					hash := hash_2(seed, x_primed, y_primed)
					idx := hash & (255 << 1)

					vec_x := f64(xi - x) + rand_vecs2d[idx] * cellular_jitter
					vec_y := f64(yi - y) + rand_vecs2d[idx | 1] * cellular_jitter

					new_distance := fast_abs(vec_x) + fast_abs(vec_y)

					distance1 = fast_max(fast_min(distance1, new_distance), distance0)
					if new_distance < distance0 {
						distance0 = new_distance
						closest_hash = hash
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
		.hybrid {
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					hash := hash_2(seed, x_primed, y_primed)
					idx := hash & (255 << 1)

					vec_x := f64(xi - x) + rand_vecs2d[idx] * cellular_jitter
					vec_y := f64(yi - y) + rand_vecs2d[idx | 1] * cellular_jitter

					new_distance := (fast_abs(vec_x) + fast_abs(vec_y)) + (vec_x * vec_x +
						vec_y * vec_y)

					distance1 = fast_max(fast_min(distance1, new_distance), distance0)
					if new_distance < distance0 {
						distance0 = new_distance
						closest_hash = hash
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
		else { // euclidean, euclidean_sq
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					hash := hash_2(seed, x_primed, y_primed)
					idx := hash & (255 << 1)

					vec_x := f64(xi - x) + rand_vecs2d[idx] * cellular_jitter
					vec_y := f64(yi - y) + rand_vecs2d[idx | 1] * cellular_jitter

					new_distance := vec_x * vec_x + vec_y * vec_y

					distance1 = fast_max(fast_min(distance1, new_distance), distance0)
					if new_distance < distance0 {
						distance0 = new_distance
						closest_hash = hash
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
	}

	if fast.m_cellular_distance_function == .euclidean
		&& int(fast.m_cellular_return_type) >= int(CellularReturnType.distance) {
		distance0 = fast_sqrt(distance0)

		if int(fast.m_cellular_return_type) >= int(CellularReturnType.distance2) {
			distance1 = fast_sqrt(distance1)
		}
	}

	match fast.m_cellular_return_type {
		.cell_value {
			return f64(closest_hash) * (1.0 / 2147483648.0)
		}
		.distance {
			return distance0 - 1
		}
		.distance2 {
			return distance1 - 1
		}
		.distance2_add {
			return (distance1 + distance0) * 0.5 - 1
		}
		.distance2_sub {
			return distance1 - distance0 - 1
		}
		.distance2_mul {
			return distance1 * distance0 * 0.5 - 1
		}
		.distance2_div {
			return distance0 / distance1 - 1
		}
	}
}

fn (fast FastNoiseLite) single_cellular_3(seed int, x f64, y f64, z f64) f64 {
	xr := fast_round(x)
	yr := fast_round(y)
	zr := fast_round(z)

	mut distance0 := 1e10 // f64
	mut distance1 := 1e10
	mut closest_hash := 0

	cellular_jitter := 0.39614353 * fast.m_cellular_jitter_modifier

	mut x_primed := int((xr - 1) * noise.prime_x)
	y_primed_base := int((yr - 1) * noise.prime_y)
	z_primed_base := int((zr - 1) * noise.prime_z)

	match fast.m_cellular_distance_function {
		.euclidean, .euclidean_sq {
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					mut z_primed := z_primed_base

					for zi := zr - 1; zi <= zr + 1; zi++ {
						hash := hash_3(seed, x_primed, y_primed, z_primed)
						idx := hash & (255 << 2)

						vec_x := f64(xi - x) + rand_vecs3d[idx] * cellular_jitter
						vec_y := f64(yi - y) + rand_vecs3d[idx | 1] * cellular_jitter
						vec_z := f64(zi - z) + rand_vecs3d[idx | 2] * cellular_jitter

						new_distance := vec_x * vec_x + vec_y * vec_y + vec_z * vec_z

						distance1 = fast_max(fast_min(distance1, new_distance), distance0)
						if new_distance < distance0 {
							distance0 = new_distance
							closest_hash = hash
						}
						z_primed += noise.prime_z
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
		.manhattan {
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					mut z_primed := z_primed_base

					for zi := zr - 1; zi <= zr + 1; zi++ {
						hash := hash_3(seed, x_primed, y_primed, z_primed)
						idx := hash & (255 << 2)

						vec_x := f64(xi - x) + rand_vecs3d[idx] * cellular_jitter
						vec_y := f64(yi - y) + rand_vecs3d[idx | 1] * cellular_jitter
						vec_z := f64(zi - z) + rand_vecs3d[idx | 2] * cellular_jitter

						new_distance := fast_abs(vec_x) + fast_abs(vec_y) + fast_abs(vec_z)

						distance1 = fast_max(fast_min(distance1, new_distance), distance0)
						if new_distance < distance0 {
							distance0 = new_distance
							closest_hash = hash
						}
						z_primed += noise.prime_z
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
		.hybrid {
			for xi := xr - 1; xi <= xr + 1; xi++ {
				mut y_primed := y_primed_base

				for yi := yr - 1; yi <= yr + 1; yi++ {
					mut z_primed := z_primed_base

					for zi := zr - 1; zi <= zr + 1; zi++ {
						hash := hash_3(seed, x_primed, y_primed, z_primed)
						idx := hash & (255 << 2)

						vec_x := f64(xi - x) + rand_vecs3d[idx] * cellular_jitter
						vec_y := f64(yi - y) + rand_vecs3d[idx | 1] * cellular_jitter
						vec_z := f64(zi - z) + rand_vecs3d[idx | 2] * cellular_jitter

						new_distance := (fast_abs(vec_x) + fast_abs(vec_y) + fast_abs(vec_z)) +
							(vec_x * vec_x + vec_y * vec_y + vec_z * vec_z)

						distance1 = fast_max(fast_min(distance1, new_distance), distance0)
						if new_distance < distance0 {
							distance0 = new_distance
							closest_hash = hash
						}
						z_primed += noise.prime_z
					}
					y_primed += noise.prime_y
				}
				x_primed += noise.prime_x
			}
		}
	}

	if fast.m_cellular_distance_function == .euclidean
		&& int(fast.m_cellular_return_type) >= int(CellularReturnType.distance) {
		distance0 = fast_sqrt(distance0)

		if int(fast.m_cellular_return_type) >= int(CellularReturnType.distance2) {
			distance1 = fast_sqrt(distance1)
		}
	}

	match fast.m_cellular_return_type {
		.cell_value {
			return f64(closest_hash) * (1 / 2147483648.0)
		}
		.distance {
			return distance0 - 1
		}
		.distance2 {
			return distance1 - 1
		}
		.distance2_add {
			return (distance1 + distance0) * 0.5 - 1
		}
		.distance2_sub {
			return distance1 - distance0 - 1
		}
		.distance2_mul {
			return distance1 * distance0 * 0.5 - 1
		}
		.distance2_div {
			return distance0 / distance1 - 1
		}
	}
}

fn (fast FastNoiseLite) single_perlin_2(seed int, x f64, y f64) f64 {
	mut x0 := fast_floor(x)
	mut y0 := fast_floor(y)

	xd0 := x - x0
	yd0 := y - y0
	xd1 := xd0 - 1
	yd1 := yd0 - 1

	xs := interp_quintic(xd0)
	ys := interp_quintic(yd0)

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y

	xf0 := lerp(grad_coord_2(seed, x0, y0, xd0, yd0), grad_coord_2(seed, x1, y0, xd1,
		yd0), xs)
	xf1 := lerp(grad_coord_2(seed, x0, y1, xd0, yd1), grad_coord_2(seed, x1, y1, xd1,
		yd1), xs)

	return lerp(xf0, xf1, ys) * 1.4247691104677813
}

fn (fast FastNoiseLite) single_perlin_3(seed int, x f64, y f64, z f64) f64 {
	mut x0 := fast_floor(x)
	mut y0 := fast_floor(y)
	mut z0 := fast_floor(z)

	xd0 := x - x0
	yd0 := y - y0
	zd0 := z - z0
	xd1 := xd0 - 1
	yd1 := yd0 - 1
	zd1 := zd0 - 1

	xs := interp_quintic(xd0)
	ys := interp_quintic(yd0)
	zs := interp_quintic(zd0)

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	z0 *= noise.prime_z
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y
	z1 := z0 + noise.prime_z

	xf00 := lerp(grad_coord_3(seed, x0, y0, z0, xd0, yd0, zd0), grad_coord_3(seed, x1,
		y0, z0, xd1, yd0, zd0), xs)
	xf10 := lerp(grad_coord_3(seed, x0, y1, z0, xd0, yd1, zd0), grad_coord_3(seed, x1,
		y1, z0, xd1, yd1, zd0), xs)
	xf01 := lerp(grad_coord_3(seed, x0, y0, z1, xd0, yd0, zd1), grad_coord_3(seed, x1,
		y0, z1, xd1, yd0, zd1), xs)
	xf11 := lerp(grad_coord_3(seed, x0, y1, z1, xd0, yd1, zd1), grad_coord_3(seed, x1,
		y1, z1, xd1, yd1, zd1), xs)

	yf0 := lerp(xf00, xf10, ys)
	yf1 := lerp(xf01, xf11, ys)

	return lerp(yf0, yf1, zs) * 0.964921414852142333984375
}

fn (fast FastNoiseLite) single_value_cubic_2(seed int, x f64, y f64) f64 {
	mut x1 := fast_floor(x)
	mut y1 := fast_floor(y)

	xs := x - x1
	ys := y - y1

	x1 *= noise.prime_x
	y1 *= noise.prime_y

	x0 := x1 - noise.prime_x
	y0 := y1 - noise.prime_y
	x2 := x1 + noise.prime_x
	y2 := y1 + noise.prime_y
	x3 := x1 + int(noise.prime_x << 1)
	y3 := y1 + int(noise.prime_y << 1)

	return cubic_lerp(cubic_lerp(val_coord_2(seed, x0, y0), val_coord_2(seed, x1, y0),
		val_coord_2(seed, x2, y0), val_coord_2(seed, x3, y0), xs), cubic_lerp(val_coord_2(seed,
		x0, y1), val_coord_2(seed, x1, y1), val_coord_2(seed, x2, y1), val_coord_2(seed,
		x3, y1), xs), cubic_lerp(val_coord_2(seed, x0, y2), val_coord_2(seed, x1, y2),
		val_coord_2(seed, x2, y2), val_coord_2(seed, x3, y2), xs), cubic_lerp(val_coord_2(seed,
		x0, y3), val_coord_2(seed, x1, y3), val_coord_2(seed, x2, y3), val_coord_2(seed,
		x3, y3), xs), ys) * (1 / (1.5 * 1.5))
}

fn (fast FastNoiseLite) single_value_cubic_3(seed int, x f64, y f64, z f64) f64 {
	mut x1 := fast_floor(x)
	mut y1 := fast_floor(y)
	mut z1 := fast_floor(z)

	xs := x - x1
	ys := y - y1
	zs := z - z1

	x1 *= noise.prime_x
	y1 *= noise.prime_y
	z1 *= noise.prime_z

	x0 := x1 - noise.prime_x
	y0 := y1 - noise.prime_y
	z0 := z1 - noise.prime_z
	x2 := x1 + noise.prime_x
	y2 := y1 + noise.prime_y
	z2 := z1 + noise.prime_z
	x3 := x1 + int(noise.prime_x << 1)
	y3 := y1 + int(noise.prime_y << 1)
	z3 := z1 + int(noise.prime_z << 1)

	return cubic_lerp(cubic_lerp(cubic_lerp(val_coord_3(seed, x0, y0, z0), val_coord_3(seed,
		x1, y0, z0), val_coord_3(seed, x2, y0, z0), val_coord_3(seed, x3, y0, z0), xs),
		cubic_lerp(val_coord_3(seed, x0, y1, z0), val_coord_3(seed, x1, y1, z0), val_coord_3(seed,
		x2, y1, z0), val_coord_3(seed, x3, y1, z0), xs), cubic_lerp(val_coord_3(seed,
		x0, y2, z0), val_coord_3(seed, x1, y2, z0), val_coord_3(seed, x2, y2, z0), val_coord_3(seed,
		x3, y2, z0), xs), cubic_lerp(val_coord_3(seed, x0, y3, z0), val_coord_3(seed,
		x1, y3, z0), val_coord_3(seed, x2, y3, z0), val_coord_3(seed, x3, y3, z0), xs),
		ys), cubic_lerp(cubic_lerp(val_coord_3(seed, x0, y0, z1), val_coord_3(seed, x1,
		y0, z1), val_coord_3(seed, x2, y0, z1), val_coord_3(seed, x3, y0, z1), xs), cubic_lerp(val_coord_3(seed,
		x0, y1, z1), val_coord_3(seed, x1, y1, z1), val_coord_3(seed, x2, y1, z1), val_coord_3(seed,
		x3, y1, z1), xs), cubic_lerp(val_coord_3(seed, x0, y2, z1), val_coord_3(seed,
		x1, y2, z1), val_coord_3(seed, x2, y2, z1), val_coord_3(seed, x3, y2, z1), xs),
		cubic_lerp(val_coord_3(seed, x0, y3, z1), val_coord_3(seed, x1, y3, z1), val_coord_3(seed,
		x2, y3, z1), val_coord_3(seed, x3, y3, z1), xs), ys), cubic_lerp(cubic_lerp(val_coord_3(seed,
		x0, y0, z2), val_coord_3(seed, x1, y0, z2), val_coord_3(seed, x2, y0, z2), val_coord_3(seed,
		x3, y0, z2), xs), cubic_lerp(val_coord_3(seed, x0, y1, z2), val_coord_3(seed,
		x1, y1, z2), val_coord_3(seed, x2, y1, z2), val_coord_3(seed, x3, y1, z2), xs),
		cubic_lerp(val_coord_3(seed, x0, y2, z2), val_coord_3(seed, x1, y2, z2), val_coord_3(seed,
		x2, y2, z2), val_coord_3(seed, x3, y2, z2), xs), cubic_lerp(val_coord_3(seed,
		x0, y3, z2), val_coord_3(seed, x1, y3, z2), val_coord_3(seed, x2, y3, z2), val_coord_3(seed,
		x3, y3, z2), xs), ys), cubic_lerp(cubic_lerp(val_coord_3(seed, x0, y0, z3), val_coord_3(seed,
		x1, y0, z3), val_coord_3(seed, x2, y0, z3), val_coord_3(seed, x3, y0, z3), xs),
		cubic_lerp(val_coord_3(seed, x0, y1, z3), val_coord_3(seed, x1, y1, z3), val_coord_3(seed,
		x2, y1, z3), val_coord_3(seed, x3, y1, z3), xs), cubic_lerp(val_coord_3(seed,
		x0, y2, z3), val_coord_3(seed, x1, y2, z3), val_coord_3(seed, x2, y2, z3), val_coord_3(seed,
		x3, y2, z3), xs), cubic_lerp(val_coord_3(seed, x0, y3, z3), val_coord_3(seed,
		x1, y3, z3), val_coord_3(seed, x2, y3, z3), val_coord_3(seed, x3, y3, z3), xs),
		ys), zs) * (1 / (1.5 * 1.5 * 1.5))
}

fn (fast FastNoiseLite) single_value_2(seed int, x f64, y f64) f64 {
	mut x0 := fast_floor(x)
	mut y0 := fast_floor(y)

	xs := interp_hermite(x - x0)
	ys := interp_hermite(y - y0)

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y

	xf0 := lerp(val_coord_2(seed, x0, y0), val_coord_2(seed, x1, y0), xs)
	xf1 := lerp(val_coord_2(seed, x0, y1), val_coord_2(seed, x1, y1), xs)

	return lerp(xf0, xf1, ys)
}

fn (fast FastNoiseLite) single_value_3(seed int, x f64, y f64, z f64) f64 {
	mut x0 := fast_floor(x)
	mut y0 := fast_floor(y)
	mut z0 := fast_floor(z)

	xs := interp_hermite(x - x0)
	ys := interp_hermite(y - y0)
	zs := interp_hermite(z - z0)

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	z0 *= noise.prime_z
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y
	z1 := z0 + noise.prime_z

	xf00 := lerp(val_coord_3(seed, x0, y0, z0), val_coord_3(seed, x1, y0, z0), xs)
	xf10 := lerp(val_coord_3(seed, x0, y1, z0), val_coord_3(seed, x1, y1, z0), xs)
	xf01 := lerp(val_coord_3(seed, x0, y0, z1), val_coord_3(seed, x1, y0, z1), xs)
	xf11 := lerp(val_coord_3(seed, x0, y1, z1), val_coord_3(seed, x1, y1, z1), xs)

	yf0 := lerp(xf00, xf10, ys)
	yf1 := lerp(xf01, xf11, ys)

	return lerp(yf0, yf1, zs)
}

// domain_warp_2 2D warps the input position using the current domain warp settings
//
// Example:
// ```
// fast.domain_warp_2(mut x, mut y)
// fast.get_noise_2(x,y)
// ```
pub fn (fast FastNoiseLite) domain_warp_2(mut x &f64, mut y &f64) {
	match fast.m_fractal_type {
		.domain_warp_progressive {
			fast.domain_warp_fractal_progressive_2(mut x, mut y)
		}
		.domain_warp_independent {
			fast.domain_warp_fractal_independent_2(mut x, mut y)
		}
		else {
			fast.domain_warp_single_2(mut x, mut y)
		}
	}
}

// domain_warp_3 3D warps the input position using the current domain warp settings
//
// Example:
// ```
// fast.domain_warp_3(mut x, mut y, mut z)
// fast.get_noise_3(x,y,z)
// ```
pub fn (fast FastNoiseLite) domain_warp_3(mut x &f64, mut y &f64, mut z &f64) {
	match fast.m_fractal_type {
		.domain_warp_progressive {
			fast.domain_warp_fractal_progressive_3(mut x, mut y, mut z)
		}
		.domain_warp_independent {
			fast.domain_warp_fractal_independent_3(mut x, mut y, mut z)
		}
		else {
			fast.domain_warp_single_3(mut x, mut y, mut z)
		}
	}
}

fn (fast FastNoiseLite) domain_warp_single_2(mut x &f64, mut y &f64) {
	seed := fast.m_seed
	amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	freq := fast.m_freqency

	mut sx := x
	mut sy := y
	fast.transform_domain_warp_coordinate_2(mut sx, mut sy)

	fast.do_single_domain_warp_2(seed, amp, freq, sx, sy, mut x, mut y)
}

fn (fast FastNoiseLite) domain_warp_single_3(mut x &f64, mut y &f64, mut z &f64) {
	seed := fast.m_seed
	amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	freq := fast.m_freqency

	mut sx := x
	mut sy := y
	mut sz := z
	fast.transform_domain_warp_coordinate_3(mut sx, mut sy, mut sz)

	fast.do_single_domain_warp_3(seed, amp, freq, sx, sy, sz, mut x, mut y, mut z)
}

fn (fast FastNoiseLite) transform_domain_warp_coordinate_2(mut x &f64, mut y &f64) {
	match fast.m_domain_warp_type {
		.open_simplex2, .open_simplex2_reduced {
			t := (*x + *y) * noise.f2
			x += t
			y += t
		}
		else {}
	}
}

fn (fast FastNoiseLite) transform_domain_warp_coordinate_3(mut x &f64, mut y &f64, mut z &f64) {
	match fast.m_warp_transform_type3d {
		.improve_xy_planes {
			xy := *x + *y
			s2 := xy * -0.211324865405187
			z *= 0.577350269189626
			x += s2 - z
			y += s2 - z
			z += xy * 0.577350269189626
		}
		.improve_zx_planes {
			xz := *x + *z
			s2 := xz * -0.211324865405187
			y *= 0.577350269189626
			x += s2 - y
			z += s2 - y
			y += xz * 0.577350269189626
		}
		.default_open_simplex2 {
			r := (*x + *y + *z) * noise.r3
			x = r - *x
			y = r - *y
			z = r - *z
		}
		else {}
	}
}

fn (fast FastNoiseLite) do_single_domain_warp_2(seed int, amp f64, freq f64, x f64, y f64, mut xr &f64, mut yr &f64) {
	match fast.m_domain_warp_type {
		.open_simplex2 {
			fast.single_domain_warp_simplex_gradient(seed, amp * 38.283687591552734375,
				freq, x, y, mut xr, mut yr, false)
		}
		.open_simplex2_reduced {
			fast.single_domain_warp_simplex_gradient(seed, amp * 16.0, freq, x, y, mut
				xr, mut yr, true)
		}
		.basic_grid {
			fast.single_domain_warp_basic_grid_2(seed, amp, freq, x, y, mut xr, mut yr)
		}
	}
}

fn (fast FastNoiseLite) do_single_domain_warp_3(seed int, amp f64, freq f64, x f64, y f64, z f64, mut xr &f64, mut yr &f64, mut zr &f64) {
	match fast.m_domain_warp_type {
		.open_simplex2 {
			fast.single_domain_warp_open_simplex2_gradient(seed, amp * 32.69428253173828125,
				freq, x, y, z, mut xr, mut yr, mut zr, false)
		}
		.open_simplex2_reduced {
			fast.single_domain_warp_open_simplex2_gradient(seed, amp * 7.71604938271605,
				freq, x, y, z, mut xr, mut yr, mut zr, true)
		}
		.basic_grid {
			fast.single_domain_warp_basic_grid_3(seed, amp, freq, x, y, z, mut xr, mut
				yr, mut zr)
		}
	}
}

fn (fast FastNoiseLite) single_domain_warp_basic_grid_2(seed int, warp_amp f64, frequency f64, x f64, y f64, mut xr &f64, mut yr &f64) {
	xf := x * frequency
	yf := y * frequency

	mut x0 := fast_floor(xf)
	mut y0 := fast_floor(yf)

	xs := interp_hermite(f64(xf - x0))
	ys := interp_hermite(f64(yf - y0))

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y

	mut hash0 := hash_2(seed, x0, y0) & (255 << 1)
	mut hash1 := hash_2(seed, x1, y0) & (255 << 1)

	lx0x := lerp(rand_vecs2d[hash0], rand_vecs2d[hash1], xs)
	ly0x := lerp(rand_vecs2d[hash0 | 1], rand_vecs2d[hash1 | 1], xs)

	hash0 = hash_2(seed, x0, y1) & (255 << 1)
	hash1 = hash_2(seed, x1, y1) & (255 << 1)

	lx1x := lerp(rand_vecs2d[hash0], rand_vecs2d[hash1], xs)
	ly1x := lerp(rand_vecs2d[hash0 | 1], rand_vecs2d[hash1 | 1], xs)

	xr += lerp(lx0x, lx1x, ys) * warp_amp
	yr += lerp(ly0x, ly1x, ys) * warp_amp
}

fn (fast FastNoiseLite) single_domain_warp_basic_grid_3(seed int, warp_amp f64, frequency f64, x f64, y f64, z f64, mut xr &f64, mut yr &f64, mut zr &f64) {
	xf := x * frequency
	yf := y * frequency
	zf := z * frequency

	mut x0 := fast_floor(xf)
	mut y0 := fast_floor(yf)
	mut z0 := fast_floor(zf)

	xs := interp_hermite(f64(xf - x0))
	ys := interp_hermite(f64(yf - y0))
	zs := interp_hermite(f64(zf - z0))

	x0 *= noise.prime_x
	y0 *= noise.prime_y
	z0 *= noise.prime_z
	x1 := x0 + noise.prime_x
	y1 := y0 + noise.prime_y
	z1 := z0 + noise.prime_z

	mut hash0 := hash_3(seed, x0, y0, z0) & (255 << 2)
	mut hash1 := hash_3(seed, x1, y0, z0) & (255 << 2)

	mut lx0x := lerp(rand_vecs3d[hash0], rand_vecs3d[hash1], xs)
	mut ly0x := lerp(rand_vecs3d[hash0 | 1], rand_vecs3d[hash1 | 1], xs)
	mut lz0x := lerp(rand_vecs3d[hash0 | 2], rand_vecs3d[hash1 | 2], xs)

	hash0 = hash_3(seed, x0, y1, z0) & (255 << 2)
	hash1 = hash_3(seed, x1, y1, z0) & (255 << 2)

	mut lx1x := lerp(rand_vecs3d[hash0], rand_vecs3d[hash1], xs)
	mut ly1x := lerp(rand_vecs3d[hash0 | 1], rand_vecs3d[hash1 | 1], xs)
	mut lz1x := lerp(rand_vecs3d[hash0 | 2], rand_vecs3d[hash1 | 2], xs)

	lx0y := lerp(lx0x, lx1x, ys)
	ly0y := lerp(ly0x, ly1x, ys)
	lz0y := lerp(lz0x, lz1x, ys)

	hash0 = hash_3(seed, x0, y0, z1) & (255 << 2)
	hash1 = hash_3(seed, x1, y0, z1) & (255 << 2)

	lx0x = lerp(rand_vecs3d[hash0], rand_vecs3d[hash1], xs)
	ly0x = lerp(rand_vecs3d[hash0 | 1], rand_vecs3d[hash1 | 1], xs)
	lz0x = lerp(rand_vecs3d[hash0 | 2], rand_vecs3d[hash1 | 2], xs)

	hash0 = hash_3(seed, x0, y1, z1) & (255 << 2)
	hash1 = hash_3(seed, x1, y1, z1) & (255 << 2)

	lx1x = lerp(rand_vecs3d[hash0], rand_vecs3d[hash1], xs)
	ly1x = lerp(rand_vecs3d[hash0 | 1], rand_vecs3d[hash1 | 1], xs)
	lz1x = lerp(rand_vecs3d[hash0 | 2], rand_vecs3d[hash1 | 2], xs)

	xr += lerp(lx0y, lerp(lx0x, lx1x, ys), zs) * warp_amp
	yr += lerp(ly0y, lerp(ly0x, ly1x, ys), zs) * warp_amp
	zr += lerp(lz0y, lerp(lz0x, lz1x, ys), zs) * warp_amp
}

fn (fast FastNoiseLite) single_domain_warp_simplex_gradient(seed int, warp_amp f64, frequency f64, xx f64, yy f64, mut xr &f64, mut yr &f64, out_grad_only bool) {
	x := xx * frequency
	y := yy * frequency

	/*
	* --- Skew moved to TransformNoiseCoordinate method ---
        * const FNfloat F2 = 0.5 * (SQRT3 - 1)
        * FNfloat s = (x + y) * F2
        * x += s; y += s
	*/

	mut i := fast_floor(x)
	mut j := fast_floor(y)
	xi := x - i
	yi := y - j

	t := (xi + yi) * noise.g2
	x0 := xi - t
	y0 := yi - t

	i *= noise.prime_x
	j *= noise.prime_y

	mut vx, mut vy := f64(0), f64(0)
	mut xo, mut yo := f64(0), f64(0)

	a := 0.5 - x0 * x0 - y0 * y0
	if a > 0 {
		aaaa := (a * a) * (a * a)
		if out_grad_only {
			grad_coord_out_2(seed, i, j, mut xo, mut yo)
		} else {
			grad_coord_dual_2(seed, i, j, x0, y0, mut xo, mut yo)
		}
		vx += aaaa * xo
		vy += aaaa * yo
	}

	c := f64(2 * (1 - 2 * noise.g2) * (1 / noise.g2 - 2)) * t +
		(f64(-2 * (1 - 2 * noise.g2) * (1 - 2 * noise.g2)) + a)
	if c > 0 {
		x2 := x0 + (2 * noise.g2 - 1)
		y2 := y0 + (2 * noise.g2 - 1)
		cccc := (c * c) * (c * c)
		if out_grad_only {
			grad_coord_out_2(seed, i + noise.prime_x, j + noise.prime_y, mut xo, mut yo)
		} else {
			grad_coord_dual_2(seed, i + noise.prime_x, j + noise.prime_y, x2, y2, mut
				xo, mut yo)
		}
		vx += cccc * xo
		vy += cccc * yo
	}

	if y0 > x0 {
		x1 := x0 + noise.g2
		y1 := y0 + (noise.g2 - 1)
		b := 0.5 - x1 * x1 - y1 * y1
		if b > 0 {
			bbbb := (b * b) * (b * b)
			if out_grad_only {
				grad_coord_out_2(seed, i, j + noise.prime_y, mut xo, mut yo)
			} else {
				grad_coord_dual_2(seed, i, j + noise.prime_y, x1, y1, mut xo, mut yo)
			}
			vx += bbbb * xo
			vy += bbbb * yo
		}
	} else {
		x1 := x0 + (noise.g2 - 1)
		y1 := y0 + noise.g2
		b := 0.5 - x1 * x1 - y1 * y1
		if b > 0 {
			bbbb := (b * b) * (b * b)
			if out_grad_only {
				grad_coord_out_2(seed, i + noise.prime_x, j, mut xo, mut yo)
			} else {
				grad_coord_dual_2(seed, i + noise.prime_x, j, x1, y1, mut xo, mut yo)
			}
			vx += bbbb * xo
			vy += bbbb * yo
		}
	}

	xr += vx * warp_amp
	yr += vy * warp_amp
}

fn (fast FastNoiseLite) single_domain_warp_open_simplex2_gradient(s int, warp_amp f64, frequency f64, xx f64, yy f64, zz f64, mut xr &f64, mut yr &f64, mut zr &f64, out_grad_only bool) {
	mut seed := s
	x := xx * frequency
	y := yy * frequency
	z := zz * frequency

	/*
	* --- Rotation moved to TransformDomainWarpCoordinate method ---
        * r3 := (2.0 / 3.0)
        * r := (x + y + z) * r3 // Rotation, not skew
        * x := r - x, y := r - y, z := r - z
	*/

	mut i := fast_round(x)
	mut j := fast_round(y)
	mut k := fast_round(z)
	mut x0 := x - i
	mut y0 := y - j
	mut z0 := z - k

	mut x_n_sign := int(-x0 - 1.0) | 1
	mut y_n_sign := int(-y0 - 1.0) | 1
	mut z_n_sign := int(-z0 - 1.0) | 1

	mut ax0 := x_n_sign * -x0
	mut ay0 := y_n_sign * -y0
	mut az0 := z_n_sign * -z0

	i *= noise.prime_x
	j *= noise.prime_y
	k *= noise.prime_z

	mut vx, mut vy, mut vz := f64(0), f64(0), f64(0)
	mut xo, mut yo, mut zo := f64(0), f64(0), f64(0)

	mut a := (0.6 - x0 * x0) - (y0 * y0 + z0 * z0)
	for l := 0; l < 2; l++ {
		if a > 0 {
			aaaa := (a * a) * (a * a)
			if out_grad_only {
				grad_coord_out_3(seed, i, j, k, mut xo, mut yo, mut zo)
			} else {
				grad_coord_dual_3(seed, i, j, k, x0, y0, z0, mut xo, mut yo, mut zo)
			}
			vx += aaaa * xo
			vy += aaaa * yo
			vz += aaaa * zo
		}

		mut b := a + 1
		mut i1 := i
		mut j1 := j
		mut k1 := k
		mut x1 := x0
		mut y1 := y0
		mut z1 := z0

		if ax0 >= ay0 && ax0 >= az0 {
			x1 += x_n_sign
			b -= x_n_sign * 2 * x1
			i1 -= x_n_sign * noise.prime_x
		} else if ay0 > ax0 && ay0 >= az0 {
			y1 += y_n_sign
			b -= y_n_sign * 2 * y1
			j1 -= y_n_sign * noise.prime_y
		} else {
			z1 += z_n_sign
			b -= z_n_sign * 2 * z1
			k1 -= z_n_sign * noise.prime_z
		}

		if b > 0 {
			bbbb := (b * b) * (b * b)
			if out_grad_only {
				grad_coord_out_3(seed, i1, j1, k1, mut xo, mut yo, mut zo)
			} else {
				grad_coord_dual_3(seed, i1, j1, k1, x1, y1, z1, mut xo, mut yo, mut zo)
			}
			vx += bbbb * xo
			vy += bbbb * yo
			vz += bbbb * zo
		}

		if l == 1 {
			break
		}

		ax0 = 0.5 - ax0
		ay0 = 0.5 - ay0
		az0 = 0.5 - az0

		x0 = x_n_sign * ax0
		y0 = y_n_sign * ay0
		z0 = z_n_sign * az0

		a += (0.75 - ax0) - (ay0 + az0)

		i += (x_n_sign >> 1) & noise.prime_x
		j += (y_n_sign >> 1) & noise.prime_y
		k += (z_n_sign >> 1) & noise.prime_z

		x_n_sign = -x_n_sign
		y_n_sign = -y_n_sign
		z_n_sign = -z_n_sign

		seed += 1293373
	}

	xr += vx * warp_amp
	yr += vy * warp_amp
	zr += vz * warp_amp
}

fn (fast FastNoiseLite) domain_warp_fractal_progressive_2(mut x &f64, mut y &f64) {
	mut seed := fast.m_seed
	mut amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	mut freq := fast.m_freqency

	for _ in 0 .. fast.m_octaves {
		mut xs := x
		mut ys := y
		fast.transform_domain_warp_coordinate_2(mut xs, mut ys)

		fast.do_single_domain_warp_2(seed, amp, freq, xs, ys, mut x, mut y)

		seed++
		amp *= fast.m_gain
		freq *= fast.m_lacunarity
	}
}

fn (fast FastNoiseLite) domain_warp_fractal_progressive_3(mut x &f64, mut y &f64, mut z &f64) {
	mut seed := fast.m_seed
	mut amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	mut freq := fast.m_freqency

	for _ in 0 .. fast.m_octaves {
		mut xs := x
		mut ys := y
		mut zs := z
		fast.transform_domain_warp_coordinate_3(mut xs, mut ys, mut zs)

		fast.do_single_domain_warp_3(seed, amp, freq, xs, ys, zs, mut x, mut y, mut z)

		seed++
		amp *= fast.m_gain
		freq *= fast.m_lacunarity
	}
}

fn (fast FastNoiseLite) domain_warp_fractal_independent_2(mut x &f64, mut y &f64) {
	mut xs := x
	mut ys := y
	fast.transform_domain_warp_coordinate_2(mut xs, mut ys)

	mut seed := fast.m_seed
	mut amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	mut freq := fast.m_freqency

	for _ in 0 .. fast.m_octaves {
		fast.do_single_domain_warp_2(seed, amp, freq, xs, ys, mut x, mut y)

		seed++
		amp *= fast.m_gain
		freq *= fast.m_lacunarity
	}
}

fn (fast FastNoiseLite) domain_warp_fractal_independent_3(mut x &f64, mut y &f64, mut z &f64) {
	mut xs := x
	mut ys := y
	mut zs := z
	fast.transform_domain_warp_coordinate_3(mut xs, mut ys, mut zs)

	mut seed := fast.m_seed
	mut amp := fast.m_domain_warp_amp * fast.m_fractal_bounding
	mut freq := fast.m_freqency

	for _ in 0 .. fast.m_octaves {
		fast.do_single_domain_warp_3(seed, amp, freq, xs, ys, zs, mut x, mut y, mut z)

		seed++
		amp *= fast.m_gain
		freq *= fast.m_lacunarity
	}
}

pub fn text_noise(c FastNoiseConfig, width int, height int, warp bool) {
	mut fast := new_noise(c)

	// println('>> $fast.m_seed')
	for xx in 0 .. height {
		for yy in 0 .. width {
			mut x, mut y := f64(xx), f64(yy)
			if warp {
				fast.domain_warp_2(mut x, mut y)
			}
			mut n := fast.get_noise_2(x, y)
			n += 1
			n /= 2
			match true {
				n < 1.0 / 17 {
					print(' ')
				}
				n < 2.0 / 17 {
					print('.')
				}
				n < 3.0 / 17 {
					print("'")
				}
				n < 4.0 / 17 {
					print(',')
				}
				n < 5.0 / 17 {
					print(';')
				}
				n < 6.0 / 17 {
					print('c')
				}
				n < 7.0 / 17 {
					print('l')
				}
				n < 8.0 / 17 {
					print('o')
				}
				n < 9.0 / 17 {
					print('d')
				}
				n < 10.0 / 17 {
					print('x')
				}
				n < 11.0 / 17 {
					print('k')
				}
				n < 12.0 / 17 {
					print('0')
				}
				n < 13.0 / 17 {
					print('O')
				}
				n < 14.0 / 17 {
					print('K')
				}
				n < 15.0 / 17 {
					print('X')
				}
				n < 16.0 / 17 {
					print('N')
				}
				else {
					print('W')
				}
			}
		}
		println('')
	}
}
