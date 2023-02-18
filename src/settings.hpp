/*
Settings for compile time options
*/

#pragma once

/*
Instantiate templates with several parameters. This is supposed to reduce
compilation time by having the required templates precompiled, but in practice
it does not help because this program uses templates everywhere.
Set to 0 to disable, 1 to enable.
*/
#define INSTANTIATE_TEMPLATES 0
