/*
Renderer for the buffer output

options:
-h,--help         show usage help message
-f,--flame        flame json file (required)
                  can use - for stdin
-i,--input        input buffers to add at start
                  (repeat option 0 or more times)
                  can use - for stdin
-o,--output       output file (required)
                  can use - for stdout
-s,--samples      number of samples to render
-t,--threads      number of threads to use for rendering
-b,--batch-size   number of samples to render per thread work unit
                  if 0 or unspecified, a reasonable size is calculated
-z,--bad-values   number of bad values allowed before terminating render
                  bad values indicate that the flame is not contractive

TODO
- support (semi) seeded random number generation (specify seed for each thread)
  - should be deterministic with 1 thread and seed specified
- support supersampling for larger in memory buffer but same size output
  - render with larger buffer but average cells into one before output
*/

#include <algorithm>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

#include <boost/program_options.hpp>

#include "renderers/buffer_renderer.hpp"

#include "utils/endian.hpp"
#include "utils/hardware.hpp"
#include "utils/json.hpp"

const std::string VERSION = "unspecified";

namespace bpo = boost::program_options;

// variables for command line arguments
size_t arg_threads = number_of_threads();
std::string arg_flame;
std::string arg_output;
std::vector<std::string> arg_input;
size_t arg_samples = 0;
size_t arg_batch_size = 0;
size_t arg_bad_values = 1 << 8;

// parsed data
tkoz::flame::Json json_flame;

// templated by dimensions to run renderer
template <size_t dims>
int run_renderer();

int main(int argc, char **argv)
{
    // setup the argument parsing with boost
    bpo::options_description options("ffr-buf usage");
    options.add_options()
        ("help,h","show this help message")
        ("flame,f",bpo::value<std::string>()->required(),"flame parameters JSON file")
        ("input,i",bpo::value<std::vector<std::string>>(),"buffers to add initially")
        ("output,o",bpo::value<std::string>()->required(),"output file")
        ("samples,s",bpo::value<size_t>()->default_value(0),"samples to render")
        ("threads,t",bpo::value<size_t>()->default_value(arg_threads),"threads to use")
        ("batch-size,b",bpo::value<size_t>()->default_value(arg_batch_size),"batch size")
        ("bad-values,z",bpo::value<size_t>()->default_value(arg_bad_values),"bad values limit");
    bpo::variables_map args_map;
    bpo::store(bpo::command_line_parser(argc,argv).options(options).run(),args_map);

    // show help message
    if (args_map.count("help") || args_map.empty() || argc < 2)
    {
        std::cerr << options;
        return 1;
    }

    // parse arguments
    arg_threads = args_map["threads"].as<size_t>();
    arg_flame = args_map["flame"].as<std::string>();
    arg_output = args_map["output"].as<std::string>();
    if (args_map.find("input") != args_map.end())
        arg_input = args_map["input"].as<std::vector<std::string>>();
    arg_samples = args_map["samples"].as<size_t>();
    arg_batch_size = args_map["batch-size"].as<size_t>();
    arg_bad_values = args_map["bad-values"].as<size_t>();

    // calculate an appropriate batch size if it is 0
    bool calculated_batch_size = false;
    if (arg_batch_size == 0)
    {
        size_t guess_batch_size = (arg_samples+255) >> 8;
        arg_batch_size = std::clamp<size_t>(guess_batch_size,1<<12,1<<20);
        calculated_batch_size = true;
    }

    // print options
    std::cerr << "ffr-buf version " << VERSION << std::endl;
    std::cerr << "--flame " << arg_flame << std::endl;
    for (auto& s : arg_input)
        std::cerr << "--input " << s << std::endl;
    std::cerr << "--output " << arg_output << std::endl;
    std::cerr << "--samples " << arg_samples;
    if (arg_samples == 0)
        std::cerr << " (not rendering)";
    std::cerr << std::endl;
    std::cerr << "--threads " << arg_threads << std::endl;
    if (calculated_batch_size)
        std::cerr << "--batch_size 0 (" << arg_batch_size << ")" << std::endl;
    else
        std::cerr << "--batch_size " << arg_batch_size << std::endl;
    std::cerr << "--bad_values " << arg_bad_values << std::endl;
    std::cerr << "--" << std::endl;

    // read and print flame data
    if (arg_flame == "-")
        json_flame = tkoz::flame::Json(std::cin);
    else
    {
        std::ifstream json_file(arg_flame,std::ios::in|std::ios::binary);
        json_flame = tkoz::flame::Json(json_file);
        json_file.close();
    }
    std::cerr << "flame: " << json_flame << std::endl;
    tkoz::flame::JsonInt dims = json_flame["dimensions"].intValue();

    // run appropriate renderer
    switch(dims)
    {
    case 1: return run_renderer<1>();
    case 2: return run_renderer<2>();
    case 3: return run_renderer<3>();
    default:
        std::cerr << "ERROR: " << dims << "D not supported" << std::endl;
        return 1;
    }
    return 0;
}

template <size_t dims>
int run_renderer()
{
    // parse flame and initialize renderer object
    tkoz::flame::Flame<dims> flame(json_flame);
    tkoz::flame::BufferRenderer<dims> renderer(flame);
    std::cerr << "buffer: " << (renderer.getBuffer().size()
        *sizeof(tkoz::flame::hist_t)) << " bytes, ";
#if TKOZ_LITTLE_ENDIAN
    std::cerr << "little endian, ";
#else
    std::cerr << "big endian, ";
#endif
    std::cerr << sizeof(tkoz::flame::hist_t) << " byte numbers" << std::endl;
    std::cerr << "size:";
    for (size_t d : flame.getSize())
        std::cerr << " " << d;
    std::cerr << std::endl;
    std::cerr << "color: " << renderer.getColorDims()
        << " dimensions" << std::endl;

    // add initial input buffers
    for (auto& s : arg_input)
    {
        std::cerr << "adding input file " << s << std::endl;
        if (s == "-")
        {
            if (!renderer.addBuffer(std::cin))
                throw std::runtime_error("error reading file");
        }
        else
        {
            std::ifstream input_file(s);
            if (!renderer.addBuffer(input_file))
                throw std::runtime_error("error reading file");
            input_file.close();
        }
    }

    // perform the render
    if (arg_samples > 0)
    {
        std::cerr << "render start" << std::endl;
        struct timespec t1,t2;
        size_t total_batches = (arg_samples + arg_batch_size - 1)
            / arg_batch_size;
        size_t batches_completed = 0;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        bool success = renderer.render(arg_samples,arg_threads,
            arg_batch_size,arg_bad_values,
            [&total_batches,&batches_completed,&t1]()
            {
                struct timespec t;
                clock_gettime(CLOCK_MONOTONIC,&t);
                size_t tn = 1000000000uLL * (t.tv_sec - t1.tv_sec)
                    + (t.tv_nsec-t1.tv_nsec);
                ++batches_completed;
                size_t tavg = tn / batches_completed;
                size_t trem = tavg * (total_batches - batches_completed);
                double p = batches_completed / (double) total_batches;
                std::cerr << "\rprogress: " << batches_completed
                    << "/" << total_batches << " batches, "
                    << std::min(arg_samples,batches_completed*arg_batch_size)
                    << "/" << arg_samples << " samples ("
                    << (int)(100*p) << "%) "
                    << (tn / 1e9) << "sec elapsed, (~"
                    << (trem / 1e9) << "sec remaining)";
            },
            [](const std::thread& t, size_t id)
            {
                std::cerr << "started thread " << id
                    << " (" << t.get_id() << ")" << std::endl;;
            });
        clock_gettime(CLOCK_MONOTONIC,&t2);
        std::cerr << std::endl;
        size_t nsecs = 1000000000uLL * (t2.tv_sec - t1.tv_sec)
            + (t2.tv_nsec - t1.tv_nsec);
        double secs = nsecs / 1e9;
        if (secs < 1e-9) // make sure it is not zero
            secs = 1e-9;
        std::cerr << "render done: " << secs << " sec ("
            << (size_t)(arg_samples/secs) << " samples/sec)" << std::endl;
        if (!success)
            std::cerr << "render failure "
                << " (flame may not be sufficiently contractive)" << std::endl;

        // output render statistics
        double iter_part = renderer.getSamplesIterated() / (double) arg_samples;
        std::cerr << "samples iterated: "
            << renderer.getSamplesIterated()
            << " (" << (100*iter_part) << "%)" << std::endl;
        double plot_part = renderer.getSamplesPlotted() / (double) arg_samples;
        std::cerr << "samples plotted: "
            << renderer.getSamplesPlotted()
            << " (" << (100*plot_part) << "%)" << std::endl;
        std::cerr << "xform selection:";
        for (size_t n : renderer.getXFormDistribution())
            std::cerr << " " << n;
        std::cerr << std::endl;
        std::cerr << "bad value xforms:";
        for (size_t n : renderer.getBadValueXForms())
            std::cerr << " " << n;
        std::cerr << std::endl;
        std::cerr << "bad value points:";
        for (auto p : renderer.getBadValuePoints())
            std::cerr << " " << p;
        std::cerr << std::endl;
        std::cerr << "extreme coordinates:";
        for (auto p : renderer.getPointExtremes())
            std::cerr << " (" << p.first << "," << p.second << ")";
        std::cerr << std::endl;
    }

    // write output
    std::cerr << "writing output" << std::endl;
    if (arg_output == "-")
    {
        if (!renderer.writeBuffer(std::cout))
            throw std::runtime_error("error writing output");
    }
    else
    {
        std::ofstream output_file(arg_output);
        if (!renderer.writeBuffer(output_file))
            throw std::runtime_error("error writing output");
        output_file.close();
    }
    return 0;
}
