# flame fractal algorithm

Flame fractals are an extension of (linear) iterated function systems. They can
be generalized to higher dimensions but to keep things simple, we consider 2D
for now. To generate a flame fractal we create some functions that map a point
in the plane to another point in the plane (more on these later).

For rendering, we divide the plane into a grid of squares and pick a random
point near $(0,0)$. Then we randomly apply functions to it and iterate a
sequence of points $x_0,x_1,x_2,\ldots$. Each time the point lands in one of
those squares on our grid, we increment a counter.

The rendering result is a probability distribution and can be plotted in
grayscale. For better detail, it should be logarithmically scaled. Color is
usually added by associating the functions with colors and measuring which grid
squares are associated with which functions.

Images can be enhanced with postprocessing such as gamma brightness adjustment
and variable radius blurring.

## functions and parameters

In $n$ dimensional space $\mathbb{R}^n$, we need some variation functions
$$ V_1,V_2,\ldots,V_m:\mathbb{R}^n\to\mathbb{R}^n $$
These functions can be defined in many different ways, even depending on
randomness. The goal is to have many options to give lots of possible variation
to flame fractals. The following are a few examples:
$$ V_1(x,y)=(x,y),\quad V_2(x,y)=(\sin(x),\cos(x)) $$
$$ V_3={1\over x^2+y^2}(x,y),\quad V_4(x,y)
={1\over\sqrt{x^2+y^2}}((x-y)(x+y),2xy) $$
Some variation functions can even have their own configurable parameters which
give even more options for variation in fractals that can be created.

Next we construct functions (also called xforms in flam3) $F_1,F_2,\ldots,F_k$
with weights $w_1,w_2,\ldots,w_k>0$ which we will assume are normalized
($w_1+w_2+\ldots+w_k=1$). For each $F_i$, we need a pre-affine function $G_i$
and post-affine function $H_i$. Affine functions are $Ax+b$ for some
$n\times n$ matrix $A$ and $n$-vector $b$. Finally, we need variation weights
$v_{i1},v_{i2},\ldots,v_{im}$. Then define
$$ F_i(x)=H_i \left( \sum_{j=1}^m v_{ij} V_j(G_i(x)) \right) $$
In practice, most $v_{ij}$ will be $0$, choosing only a small number of
variation functions to use.

Optionally, we can also define a final xform $F_f$ similarly, but it does not
have an associated weight.

For color, we assign $r$-dimensional color vectors to the xforms:
$$ c_1,c_2,\ldots,c_k\in[0,1]^r $$
There can also optionally be a final color $c_f$. We also choose a color speed
parameter $0\leq s\leq1$ which is used for exponential averaging. The color
vectors can be mapped to a color space in many different ways. A few are

- map 1D values to a color palette
- map 2D values to hue and saturation
- map 3D values to red, green, and blue

Usually the brightness is determined from the logarithmically scaled probability
distribution, but there can be other ways to use that information.

Full xform weight information can be found by using $r=k$ and choosing each
color vector to be a standard basis vector ($c_i=e_i$). The color iteration is
$$ c = sc_i + (1-s)c $$
which converges to a vector whose component sum (or 1-norm) is $1$.

The many parameters allow a lot of freedom in defining flame fractals. In
practice, we need the functions to be contractive on average to generate a good
image. Poor choice of parameters will result in degenerate images. Finding good
images is a matter of running Monte-Carlo simulations for various parameters and
selecting those that that have aesthetic value. One way to easily identify some
bad images is when the point iterates outside a large chosen bound often.

## rendering algorithm

To render a flame fractal, we first zero-initialize a buffer corresponding to a
grid in $\mathbb{R}^n$. For each cell in the grid (which would be pixels in 2D),
we have a counter ($\alpha$) and a color vector ($\beta$).

1. Choose a random point $x\in[-1,1]^n$, the unit cube.
2. Choose a random color $c\in[0,1]^r$.
3. Repeat the following loop for sufficiently many plotted samples.
    1. Choose an xform $F_i$ (where $1\leq i\leq k$).
    2. Set $x$ to $F_i(x)$.
    3. Set $c$ to $sc_i+(1-s)c$ if $F_i$ has a color $c_i$.
    3. Set $x'=F_f(x)$ if a final xform is present, otherwise $x'=x$.
    4. Set $c'=sc_f+(1-s)c$ if a final color is present, otherwise $c'=c$.
    4. Plot $x'$ and $c'$ except during the first few iterations.
        - Increment the counter ($\alpha$) for the cell $x'$ lands in.
        - Add the color vector $c'$ to the cell's color value ($\beta$).

Skipping plotting in the first few iterations allows the generated sequence to
converge to the actual solution of the iterated function system.

The colors are blended with exponential averaging which mixes the colors
according to the most recently used xforms. This allows cells to be colored
based on how much the xform contributes to the likelihood of the point landing
in that cell.

To save the full color blending information as a combination of each variation
color, we can use $c_1=e_1,c_2=e_2,\ldots,c_k=e_k$ where
$e_1,e_2,\ldots,e_k\in\mathbb{R}^k$ form the standard basis. This does require
choosing a color speed parameter. Being able to generalize from the color speed
requires using a lot more memory. Flames generally have a small number of xforms
so using $k$ dimensions for color information is not too bad.

## creating a 2d image

Here, we describe some ways to create an image from a rendering of a 2D flame
fractal. Each cell in the 2D buffer would correspond to a single pixel, or they
can be downsampled to combine multiple rendering cells into a single pixel for
improved quality. There is a lot of complexity in what is possible with creating
images from the rendering buffer.

The methods described here may differ from that of the original flam3 software.

Usually brightness is determined from the counter. Because of the high range,
the most aesthetic value is brought out by log scaling the counter to show
differences in density at various scales, which brings out more detail in an
image. The densest parts are usually exponentially more dense than the lower
density parts. Some other scaling function ideas are linear and nth roots, or
even binary.

Let a pixel have counter $\alpha$ and color vector $\beta$ (the color vector is
the average color of iterations that reached that pixel). Let $\alpha'$ be the
maximum $\alpha$ over all the pixels. Then we scale the $\alpha$ value as
$b=\log(1+\alpha)/\log(1+\alpha')\in[0,1]$ for the brightness parameter.
Next, let $b'=b^{1/\gamma}\in[0,1]$ be the brightness adjusted with a gamma
parameter where $\gamma>0$. Now, there are several methods to determine color.
A few possibilities are:

1. If $\beta$ is 1D, use a color palette to select a color and determine
brightness with $b'$.
2. If $\beta$ is 2D use the HSL/HSV color spaces. Determine hue and saturation
from $\beta$ and lightness/value from $b'$.
3. If $\beta$ is 3D, use $\beta$ as an RGB color and multiply it by $b'$ to
adjust brightness.

Lower density parts of the image can also benefit from blurring, by using a blur
radius that varies depending on how dense part of the buffer is. Another option
is oversampling and determining density of pixels as an average of several
histogram points. The lower density parts of the fractal get smoother as the
number of rendering samples increases.

## compile time constants

Real number type `tkoz::flame::num_t`: either `float` or `double`. The type used
for point iteration and calculations.

Histogram counter type `tkoz::flame::hist_t`: either `uint32_t` or `uint64_t`.
The type used for histogram frequency counting. This must be the same size in
bytes as `num_t` due to memory alignment design of the rendering buffer.

Small epsilon to avoid division by zero `tkoz::flame::eps`: $10^{-10}$ for
`float` and $10^{-20}$ for `double`. flam3 uses `1e-10`.

Machine epsilon `tkoz::flame::emach`: $2^{-23}$ for `float` and $2^{-52}$ for
`double`. This is the distance between 1 and the smallest representable number
greater than 1.

Maximum dimension `tkoz::flame::max_dim`: 65535, limits the length of the buffer
in pixels/buckets along any dimension. This limit may make sense to change.

Maximum coordinate of rendering box `tkoz::flame::max_rect`: absolute values of
rendering bounds must be less than this limit. This may make sense to change.
Currentnly $10^5$ for `float` and $10^{10}$ for `double`.

Settle iterations `tkoz::flame::settle_iters`: number of iterations to skip for
iteration to converge to the attractor. 24 for `float` and 53 for `double`.
flam3 uses 20. This value make make sense to change but should not be too big.
It may make sense to have this value depend on the number of xforms or
dimensions.

Bad value threshold `tkoz::flame::bad_value_threshold`: threshold for which a
point is assumed to diverge to infinity. $10^{10}$ for `float` and $10^{20}$
for `double`. flam3 uses `1e10`.

## sources

Original paper from Scott Draves:
https://flam3.com/flame_draves.pdf
