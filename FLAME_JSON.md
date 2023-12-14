# flame json format

Flame fractals are specified with JSON files with comments allowed.

## flame object (root)

An object with the following keys:
- `"dimensions"`: an integer $n\in[1,10]$. The C++ classes have to be compiled
for each value used. 10 is an arbitrary limit chosen for practical reasons.
- `"size"`: list of $n$ positive integers $[N_1,N_2,\ldots,N_n]$ for the number
of pixels/buckets in each dimension. Each $N_i$ must be below $2^{16}$. This
limit may make sense to increase for 2D or 3D. The product of all must be below
$2^{64}$, although 64 bit address space means the limit is much smaller
(no one has enough memory to get close to this limit anyway).
- `"bounds"`: list of $n$ lists $[B_{i,0},B_{i,1}]$ specifying the range for
plotting points on dimension $i$ as $B_{i,0}\leq x_i\leq B_{i,1}$
(where $x\in\mathbb{R}^n$)
- `"color_dimensions"`: size of color vectors. In practice, this might be 1
(single dimension color palette) or 2 (hue+sat for HSL color space) or 3
(full color per xform), but there is also reason to make this equal to the
number of xforms. These are only some ideas, there may be more useful ones.
If not specified or is 0, color is disabled. An arbitrary limit of 127 is chosen
for this parameter.
- `"color_speed"`: a paramater $s\in[0,1]$ used as the exponential averaging
parameter for blending colors. 0 leaves the initial random color unchanged and
1 means to always use the color of the current xform. If not specified, the
default used is 0.5.
- `"xforms"`: a list of xform objects (described below), must have at least 1.
- `"final_xform"`: an xform object with a few differences. Including a final
xform is optional.

## xform object

An object with the following keys:
- `"weight"`: a floating point number describing the relative probability of
xform selection. 0 means to ignore the xform (useful for testing).
- `"color"`: color vector with each component in $[0,1]$. If not specified, this
xform will not change the iterating color.
- `"variations"`: a list of variation objects (described below).
- `"pre_affine"`: an object `{"A":[...],"b":[...]}` where `"A"` is a $n\times n$
array of floating point numbers and $b$ is a $n$-vector of floating point
numbers. These describe the affine transformation $Ax+b$.
- `"post_affine"`: same format as `"pre_affine"`

## variation object

Every variation has a `"name"` which is a string for the name of the variation
and a `"weight"` which is a floating point number. Some variations have
additional parameters described below.
