/*
Renderer for buffers of variable dimension size.
The upper limit of dimensions is specified at compile time.
*/

#include <fstream>
#include <iostream>
#include <string>

#include <boost/program_options.hpp>

#include "renderers/helpers.hpp"
#include "renderers/histogram_renderer.hpp"
#include "types/types.hpp"
#include "utils/hardware.hpp"
#include "utils/misc.hpp"

const std::string VERSION = "unspecified";

namespace bpo = boost::program_options;
typedef float num_t; // can also be double (slower)
typedef uint32_t hist_t; // can also be uint64_t

int main(int argc, char **argv)
{
    size_t num_threads = number_of_threads();
    bpo::options_description options("ffr-buffer usage");
    options.add_options()
        ("help,h","show this message")
        ("flame,f",bpo::value<std::string>()->required(),
            "flame parameters JSON file (required)")
        ("output,o",bpo::value<std::string>()->required(),
            "output file (required)")
        ("input,i",bpo::value<std::vector<std::string>>(),
            "buffers to add for initial histogram, "
            "can use multiple times (default none)")
        ("samples,s",bpo::value<size_t>()->default_value(0),
            "samples to render (default 0)")
        ("threads,T",bpo::value<size_t>()->default_value(num_threads),
            "number of threads to use (default number of threads available)")
        ("batch_size,z",bpo::value<size_t>()->default_value(1<<18),
            "multithreading batch size (default 2^18)")
        ("bad_values,B",bpo::value<size_t>()->default_value(10),
            "bad value limit for terminating render (default 10)");
    bpo::variables_map args;
    bpo::store(bpo::command_line_parser(argc,argv).options(options).run(),args);
    if (args.count("help") || args.empty() || argc < 2)
    {
        std::cerr << options;
        return 1;
    }
    if (!args.count("flame")) // required argument
    {
        std::cerr << "ERROR: flame JSON not specified" << std::endl;
        return 1;
    }
    if (!args.count("output")) // required argument
    {
        std::cerr << "ERROR: output file not specified" << std::endl;
        return 1;
    }
    // cmdline values
    std::string arg_flame = args["flame"].as<std::string>();
    std::vector<std::string> args_input;
    if (args.count("input") > 0)
        args_input = args["input"].as<std::vector<std::string>>();
    size_t arg_samples = args["samples"].as<size_t>();
    std::string arg_output = args["output"].as<std::string>();
    size_t arg_threads = args["threads"].as<size_t>();
    size_t arg_batch_size = args["batch_size"].as<size_t>();
    size_t arg_bad_values = args["bad_values"].as<size_t>();
    if (arg_threads < 1 || arg_threads > 256)
    {
        std::cerr << "ERROR: must use between 1 and 256 threads" << std::endl;
        return 1;
    }
    if (arg_batch_size < (1<<8) || arg_batch_size > (1<<30))
    {
        std::cerr << "ERROR: batch size < 2^8 or > 2^30" << std::endl;
        return 1;
    }
    if (arg_bad_values < 1)
    {
        std::cerr << "ERROR: bad values must be positive" << std::endl;
        return 1;
    }
    std::cerr << "ffr-buffer version " << VERSION << std::endl;
    #if TKOZ_LITTLE_ENDIAN
    std::cerr << "endianness: little" << std::endl;
#else
    std::cerr << "endianness: big" << std::endl;
#endif
    std::cerr << "number of cpus: " << number_of_threads() << std::endl;
    // print options
    std::cerr << "command line arguments:" << std::endl;
    std::cerr << "ffbuf" << std::endl;
    std::cerr << "--flame: " << arg_flame << std::endl;
    for (std::string arg_input : args_input)
        std::cerr << "--input: " << arg_input << std::endl;
    std::cerr << "--samples: " << arg_samples << std::endl;
    std::cerr << "--output: " << arg_output << std::endl;
    std::cerr << "--threads: " << arg_threads << std::endl;
    std::cerr << "--batch_size: " << arg_batch_size << std::endl;
    std::cerr << "--bad_values: " << arg_bad_values << std::endl;
    std::cerr << "--" << std::endl;
    if (arg_threads > num_threads)
        std::cerr << "WARN: using more threads than the system supports"
            << std::endl;
    Json json_flame;
    if (arg_flame == "-")
        json_flame = Json(std::cin);
    else
        json_flame = Json(read_text_file(arg_flame));
    std::cerr << "flame json (comments removed): " << json_flame << std::endl;
    tkoz::flame::HistogramRendererInterface<num_t,hist_t> *renderer =
        tkoz::flame::instantiateHistogramRenderer<num_t,hist_t>(json_flame);
    if (!renderer)
        throw std::runtime_error("renderer instantiation failed");
    hist_t *buf = nullptr;
    for (std::string arg_input : args_input)
    {
        if (!buf)
            buf = new hist_t[renderer->getHistogramSize()];
        std::ifstream ifile;
        if (arg_input != "-")
            ifile = std::ifstream(arg_input,std::ios::in|std::ios::binary);
        std::istream& istream = arg_input == "-" ? std::cin : ifile;
        if (!istream.read((char*)buf,renderer->getHistogramSizeBytes()))
            throw std::runtime_error("failed reading "+arg_input);
        if (istream.peek() != EOF)
            std::cerr << "WARN: did not reach EOF of " << arg_input
                << std::endl;
        renderer->addHistogram(buf);
    }
    if (buf)
        delete[] buf;
    if (arg_samples)
    {
        struct timespec t1,t2;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        i32 prev_percent = -1;
        size_t prev_tsec = t1.tv_sec;
        renderer->render(arg_samples,arg_threads,arg_batch_size,arg_bad_values,
            [&prev_percent,&prev_tsec,&t1,&t2](float p)
            {
                clock_gettime(CLOCK_MONOTONIC,&t2);
                size_t nsecs = 1000000000uLL*(t2.tv_sec-t1.tv_sec)
                    +(t2.tv_nsec-t1.tv_nsec);
                size_t secs = nsecs / 1000000000;
                i32 percent = (i32)(100.0*p);
                size_t tsec = t2.tv_sec;
                if (percent > prev_percent || tsec > prev_tsec)
                {
                    std::cerr << '\r' << "rendering... " << percent << "% ("
                        << secs << " sec elapsed)";
                    prev_percent = percent;
                    prev_tsec = tsec;
                }
            });
        std::cerr << '\r' << "rendering... 100%" << std::endl;
        clock_gettime(CLOCK_MONOTONIC,&t2);
        size_t nsecs = 1000000000uLL*(t2.tv_sec-t1.tv_sec)
            +(t2.tv_nsec-t1.tv_nsec);
        float secs = (float)nsecs / 1000000000.0;
        std::cerr << "render done" << std::endl;
        fprintf(stderr,"time (seconds): %f\n",secs);
        fprintf(stderr,"samples iterated: %lu\n",
            renderer->getSamplesIterated());
        fprintf(stderr,"samples plotted: %lu\n",
            renderer->getSamplesPlotted());
        fprintf(stderr,"samples/sec: %lf\n",
            renderer->getSamplesIterated()/secs);
        double ratio = (double)renderer->getSamplesPlotted()
            / renderer->getSamplesIterated();
        fprintf(stderr,"plotted/iterated: %lf%%\n",100*ratio);
    }
    std::ofstream ofile;
    if (arg_output != "-")
        ofile = std::ofstream(arg_output,std::ios::out|std::ios::binary);
    std::ostream& ostream = arg_output == "-" ? std::cout : ofile;
    if (!renderer->writeHistogram(ostream))
        throw std::runtime_error("failed to write output");
    if (arg_output != "-")
        ofile.close();
    return 0;
}
