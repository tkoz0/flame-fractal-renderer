# flame fractal renderer

## what are flame fractals?

Flame fractals are fractals created by certain iterated function systems. They
were created by Scott Draves in 1992. See his paper here:
https://flam3.com/flame_draves.pdf. He also has a software called flam3 which is
the main inspiration for writing this renderer. You can find it here:
https://github.com/scottdraves/flam3.

Flame fractals are based on iterated function systems. They are typically
considered in 2 dimensions for an image, but can be generalized to higher
dimensions. In $n$ dimensions, we have a set of functions $f_i:R^n\to R^n$ with
weights $w_i$ that sum to 1. The details in specifying these functions will come
later. To render a flame fractal, start with a random point $x$ in the biunit
square/cube/hypercube centered at the origin. Then to perform an iteration,
randomly choose a function $f_i$ according to the weights and update $x$ to be
$f_i(x)$. This will generate a sequence of points. For each point (except some
at the beginning to allow the sequence to settle), assign it to an area/volume
in the space $R^n$ divided into tiny squares/cubes/hypercubes.

In 2D, we have functions $f_i:R^2\to R^2$. We divide the 2D plane into small
squares corresponding to pixels. For each pixel, we have a counter that gets
incremented each time the generated sequence lands inside its area. Finally, to
plot the result, the value for each pixel is scaled and then output as a
grayscale image. Currently, my renderer supports 3 methods: binary, linear,
logarithmic (default).

Flame fractals can be rendered with color. The common ways to do this are by
assigning colors to each function and blending the colors of the most recently
selected functions using exponential weighting. My renderer currently does not
support color.

The functions $f_i$ are defined by a pre-affine and post-affine transformation.
Each of these is a linear transformation that can be defined as $Ax+b$ for some
$n\times n$ matrix $A$ and $n$-vector $b$. Between these, we apply a sum of
variation functions. Then to apply $f_i$ to a point $x$, we first apply the
pre-affine transformation. Next, use this point to compute the variation
function sum. Finally, apply the post-affine transformation to this sum. 

## why another renderer when others exist?

I initially created this because I wanted to create the probability distribution
buffer and then render images based on that (since the probability distribution
buffer is the most computationally expensive part to generate, it makes sense to
save it and then do the image conversion/filtering later).

It expanded a bit when I decided to write a whole renderer that I would try to
optimize decently, although that is difficult to do when trying to support
arbitrary functions to be specified in JSON files. I ended up switching from C
to C++ and making the whole thing templated to be able to support different
precision (currently single precision is sufficient), pixel counter bit size
(there is not much reason to go past 32 bit), and random number generator (I
decided to go with ISAAC as flam3 uses ISAAC, but the code for a
java.util.Random implementation is still included).

Another reason was interest in exploring 3D flame fractals. This software has
not been extended to support 3D yet. Currently the code can support higher
dimensions but there is not yet a renderer that can make use of it.

## json format

Flames fractals for this software are specified in JSON. See the example flames
for a good idea of how to write them. This documentation should be expanded
eventually. The software does not give good error messages for incorrect format.

## building

Simply run `make` in the directory. It will build the executable `bin/ffbuf`.
This was only tested on Ubuntu 20.04 and 22.04.

## command line options

The help message:

```
ffbuf: usage
[-h --help]: show this message
-f --flame: flame parameters JSON file (required)
-o --output: output file (required)
[-i --input]: buffer (default none, render new buffer)
[-s --samples]: samples to render (default 0)
[-t --type]: output type (png,pgm,buf) (default use file extension)
[-b --img_bits]: bit depth for png or pgm output (8 or 16) (default 8)
[-T --threads]: number of threads to use (default 1)
[-z --batch_size]: multithreading batch size (default 250000)
[-B --bad_values]: bad value limit for terminating render (default 10)
[-m --scaler]: scaling function for image (binary,linear,log) (default log)
```

You can use `-` for stdin/stdout with `-f`, `-o`, and `-i`, but should not use
it for both `-f` and `-i` in the same command.

In more detail:

`-f --flame` specifies the JSON file with the flame fractal parameters

`-o --output` specifies the output file to write. The type of output is
determined automatically from the extension (.png, .pgm, or .buf, see the `-t`
option)

`-i --input` specify a buffer file, otherwise start with an empty (zeroed)
buffer. This can be used to add samples to an existing buffer or convert an
existing buffer into an image.

`-s --samples` number of samples to render

`-t --type` type of output (png, pgm, or buf)

`-b --img_bits` bit depth of image output (8 or 16)

`-T --threads` number of threads for multithreading

`-z ---batch_size` number of samples to render per work unit on a thread, higher
means less multithreading overhead

`-B --bad_values` number of bad values before considering a render to be bad
(bad values are when the points go past a certain bound that suggests they are
growing toward infinity)

`-m --scaler` scaling function for image (binary, linear, or log). Logarithmic
is recommended by Draves's paper

## example commands

Render a flame with 1000000 samples and save as an image

`ffbuf -f flame.json -o flame.png -s 1000000`

Save the probability distribution buffer instead of an image

`ffbuf -f flame.json -o flame.buf -s 1000000`

Add additional rendered samples to an existing buffer

`ffbuf -f flame.json -i flame1.buf -o flame2.buf -s 1000000`

Convert the probability distribution buffer to an image

`ffbuf -f flame.json -i flame.buf -o flame.png`

Multithreaded render with linear scaling instead of logarithmic

`ffbuf -f flame.json -o flame.png -s 1000000 -T 2 -m linear`

## explanation of the output

Rendering statistics are written to standard error. This section needs to be
expanded.
