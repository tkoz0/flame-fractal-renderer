# flame fractal algorithm

## defining the functions

Flame fractals are an extension of (linear) iterated function systems. Each
flame fractal is defined by a number of dimensions $n$ and $k$ xforms, or
transformation functions $F_1,F_2,\ldots,F_k:\mathbb{R}^n\to\mathbb{R}^n$, we
require $k\geq1$. Optionally, there may be a final xform $F_0$. First we
describe how to create these xforms $F_0,F_1,F_2,\ldots,F_k$.

We first have a set of (possibly nonlinear) variation functions. There may be
arbitrarily many and we can always create more. We will represent them as
$V_1,V_2,\ldots,V_m:\mathbb{R}^n\to\mathbb{R}^n$ supposing we have created $m$
variation functions. A few examples in 2D are:
$$V_0(x,y)=(x,y)$$
$$V_1(x,y)=(\sin(x),\sin(y))$$
$$V_2(x,y)=\frac{1}{r^2}(x,y),\quad r=\sqrt{x^2+y^2}$$

Now, for each xform $F_i$, we create them as a composition of a pre-affine
function, followed by a sum of variation functions, followed by a post-affine
function. Let $g_i$ and $h_i$ be affine functions, that is, expressible as
$Ax+b$ for some $n\times n$ matrix $A$ and $n$-vector $b$. Then an xform is
defined as
$$F_i(x)=h_i\left(\sum_{j=1}^{m}v_{ij}V_j(g_i(x))\right)$$
where each $v_{ij}\in\mathbb{R}$ are weights chosen as a representation of how
much to use the variation function. In practice, many will be $0$, and a small
number of variations will be selected with nonzero weights to be used. Each
$g_i$ and $h_i$ can be arbitrary affine functions.

Additionally, we assign probabilities $w_i\geq0$ to the xforms
$F_1,F_2,\ldots,F_k$. These describe relative weights for selecting xforms and
are normalized to sum to $1$.

For color, we assign color vectors $c_1,c_2,\ldots,c_k$ to each of the xforms,
and possibly a final color $c_0$ which may be specified whether or not there is
a final xform. The original flame fractal algorithm uses 1 dimension for these,
but we allow this to vary for more possibilities. Each coordinate on these
vectors should be in $[0,1]$. Let $r$ be the number of dimensions in the color
vectors. Also let $s$ be a color speed parameter for exponential averaging where
$s\in[0,1]$. To generalize, we allow color specification to be optional and for
each xform to have its own color speed parameter.

The many parameters allow a lot of freedom in defining flame fractals. In
practice, we need the functions to be contractive on average to generate a good
image. Poor choice of parameters will result in degenerate images. Finding good
images is a matter of running Monte-Carlo simulations for various parameters and
selecting those that that have aesthetic value.

## the rendering algorithm

To render a flame fractal, we first zero-initialize a buffer corresponding to a
grid in $\mathbb{R}^n$. For each cell in the grid (which would be pixels in 2D),
we have a counter ($\alpha$) and a color vector ($\beta$).

1. Choose a random point $x\in[-1,1]^n$, the unit cube.
2. Choose a random color $c\in[0,1]^r$.
3. Repeat the following loop for sufficiently many plotted samples.
    1. Choose an xform $F_i$ (where $1\leq i\leq k$).
    2. Set $x$ to $F_i(x)$.
    3. Set $c$ to $sc_i+(1-s)c$ if $F_i$ has a color $c_i$.
    3. Set $x'=F_0(x)$ if a final xform is present, otherwise $x'=x$.
    4. Set $c'=sc_0+(1-s)c$ if a final color is present, otherwise $c'=c$.
    4. Plot $x'$ and $c'$ except during the first few iterations.
        - Increment the counter ($\alpha$) for the cell $x'$ lands in.
        - Add the color vector $c'$ to the cell's color value.

Skipping plotting in the first few iterations allows the generated sequence to
approach the actual solution of the iterated function system.

The colors are blended with exponential averaging which mixes the colors
according to the most recently used xforms. This allows cells to be colored
based on how much the xform contributes to the likelihood of the point landing
in that cell.

To save the full color blending information as a combination of each variation
color, we can use $c_1=e_1,c_2=e_2,\ldots,c_k=e_k$ where
$e_1,e_2,\ldots,e_k\in\mathbb{R}^k$ form the standard basis. This does require
choosing a color speed parameter. Being able to generalize from the color speed
requires using a lot more memory. Flames generally have a small number of xforms
so using $k$ dimensions for color information is not bad.

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
density parts. Some other scaling function ideas are linear and nth roots.

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
