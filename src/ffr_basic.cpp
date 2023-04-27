/*
Renderer for 2d buffers and 2d images.

ffr-basic: usage
[-h --help]: show this message
-f --flame: flame parameters JSON file (required)
-o --output: output file (required)
[-i --input]: buffers to start with (can use multiple input files)
[-s --samples]: samples to render (default 0)
[-t --type]: output type (png,pgm,buf) (default use file extension)
[-b --img_bits]: bit depth for png or pgm output (8 or 16) (default 8)
[-T --threads]: number of threads to use (default 1)
[-z --batch_size]: multithreading batch size (default 250000)
[-B --bad_values]: bad value limit for terminating render (default 10)
[-m --scaler]: scaling function for image (binary,linear,log) (default log)

planned options (not available yet):
[-p --precision]: calculation precision (single or double) (default single)
[-w --hist_bits]: histogram integer bit size (32 or 64) (default 32)
[-R --rng]: random number generator (java,isaac32,isaac64) (default isaac32)
[-r --seed]: random number generator seed (default random)
[-v --verbose]: show extra information
*/

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>

#include "renderers/histogram_renderer.hpp"
#include "utils/image.hpp"
#include "utils/misc.hpp"
#include "utils/endian.hpp"
#include "utils/hardware.hpp"

const std::string VERSION = "unspecified";

namespace bpo = boost::program_options;
typedef float num_t; // can also be double (slower)
// half precision may be ok for smaller buffers
typedef uint32_t hist_t; // can also be uint64_t
// uint16_t is probably too small for practical use cases

int main(int argc, char **argv)
{
    size_t num_threads = number_of_threads();
    bpo::options_description options("ffr-basic usage");
    options.add_options()
        ("help,h","show this message")
        ("flame,f",bpo::value<std::string>()->required(),
            "flame parameters JSON file (required)")
        ("output,o",bpo::value<std::string>()->required(),
            "output file (required)")
        ("input,i",bpo::value<std::vector<std::string>>(),
            "buffers to add for initial histogram (default none)")
        ("samples,s",bpo::value<size_t>()->default_value(0),
            "samples to render (default 0)")
        ("type,t",bpo::value<std::string>()->default_value(""),
            "output type (png,pgm,buf) (default use file extension)")
        ("img_bits,b",bpo::value<size_t>()->default_value(8),
            "bit depth for png or pgm output (8 or 16) (default 8)")
        ("threads,T",bpo::value<size_t>()->default_value(num_threads),
            "number of threads to use (default number of threads available)")
        ("batch_size,z",bpo::value<size_t>()->default_value(1<<18),
            "multithreading batch size (default 2^18)")
        ("bad_values,B",bpo::value<size_t>()->default_value(10),
            "bad value limit for terminating render (default 10)")
        ("scaler,m",bpo::value<std::string>()->default_value("log"),
            "scaling function for image render "
            "(binary,linear,log) (default log)");
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
    std::string arg_type = args["type"].as<std::string>();
    size_t arg_img_bits = args["img_bits"].as<size_t>();
    size_t arg_threads = args["threads"].as<size_t>();
    size_t arg_batch_size = args["batch_size"].as<size_t>();
    size_t arg_bad_values = args["bad_values"].as<size_t>();
    std::string arg_scaler = args["scaler"].as<std::string>();
    // check args
    if (arg_type != "" && arg_type != "png" && arg_type != "pgm"
        && arg_type != "buf")
    {
        std::cerr << "ERROR: invalid output type" << std::endl;
        return 1;
    }
    if (arg_img_bits != 8 && arg_img_bits != 16)
    {
        std::cerr << "ERROR: invalid image bits per pixel" << std::endl;
        return 1;
    }
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
    if (arg_scaler != "binary" && arg_scaler != "linear" && arg_scaler != "log")
    {
        std::cerr << "ERROR: scaler must be binary/linear/log" << std::endl;
        return 1;
    }
    if (arg_type == "") // find type from extension
    {
        if (string_ends_with(arg_output,".png"))
            arg_type = "png";
        else if (string_ends_with(arg_output,".pgm"))
            arg_type = "pgm";
        else if (string_ends_with(arg_output,".buf"))
            arg_type = "buf";
        else
        {
            std::cerr << "ERROR: unable to infer output type" << std::endl;
            return 1;
        }
    }
    std::cerr << "ffr-basic version " << VERSION << std::endl;
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
    std::cerr << "--type: " << arg_type << std::endl;
    std::cerr << "--img_bits: " << arg_img_bits << std::endl;
    std::cerr << "--threads: " << arg_threads << std::endl;
    std::cerr << "--batch_size: " << arg_batch_size << std::endl;
    std::cerr << "--bad_values: " << arg_bad_values << std::endl;
    std::cerr << "--" << std::endl;
    if (arg_threads > num_threads)
        std::cerr << "WARN: using more threads than the system supports"
            << std::endl;
    // parse flame file
    Json json_flame;
    if (arg_flame == "-") // input from stdin
        json_flame = Json(std::cin);
    else
        json_flame = Json(read_text_file(arg_flame));
    std::cerr << "flame json (comments removed): " << json_flame << std::endl;
    tkoz::flame::HistogramRenderer<num_t,2,hist_t,false> renderer(json_flame);
    const tkoz::flame::Flame<num_t,2>& flame = renderer.getFlame();
    std::cerr << "x size: " << flame.getSize()[0] << std::endl;
    std::cerr << "y size: " << flame.getSize()[1] << std::endl;
    std::pair<num_t,num_t> xb = flame.getBounds()[0];
    std::pair<num_t,num_t> yb = flame.getBounds()[1];
    num_t xdiff = xb.second - xb.first;
    num_t ydiff = yb.second - yb.first;
    fprintf(stderr,"rect ratio (render bounds): %f\n",(float)ydiff/xdiff);
    fprintf(stderr,"size ratio (buffer): %f\n",
        (float)flame.getSize()[1]/flame.getSize()[0]);
    hist_t *fbuf = nullptr; // buffer for storing file contents
    // load buffer if specified
    size_t i = 0;
    for (std::string arg_input : args_input)
    {
        ++i;
        fprintf(stderr,"reading input file %s (%lu/%lu)\n",
            arg_input.c_str(),i,args_input.size());
        if (!fbuf) // use buffer to read from files
            fbuf = new hist_t[renderer.getHistogramSize()];
        if (arg_input == "-") // input from stdin
        {
            if (!std::cin.read((char*)fbuf,renderer.getHistogramSizeBytes()))
            {
                std::cerr << "ERROR: unable to read (enough) buffer bytes"
                    << std::endl;
                return 1;
            }
            if (std::cin.peek() != EOF)
                std::cerr << "WARN: did not reach EOF of stdin" << std::endl;
        }
        else
        {
            std::ifstream ifs(arg_input,std::ios::in|std::ios::binary);
            if (!ifs.read((char*)fbuf,renderer.getHistogramSizeBytes()))
            {
                std::cerr << "ERROR: unable to read (enough) buffer bytes"
                    << std::endl;
                return 1;
            }
            if (ifs.peek() != EOF)
                std::cerr << "WARN: did not reach EOF of file" << std::endl;
            ifs.close();
        }
        fprintf(stderr,"read %lu bytes from %s\n",
            renderer.getHistogramSizeBytes(),arg_input.c_str());
        renderer.addHistogram(fbuf);
    }
    if (fbuf) // clean up file buffer
        delete[] fbuf;
    size_t buffer_sum_initial = renderer.histogramSum();
    fprintf(stderr,"buffer sum (initial): %lu\n",buffer_sum_initial);
    if (arg_samples)
    {
        std::cerr << "render start" << std::endl;
        struct timespec t1,t2;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        i32 prev_percent = -1;
        size_t prev_tsec = t1.tv_sec;
        renderer.render(arg_samples,arg_threads,
            arg_batch_size,arg_bad_values,
            [&prev_percent,&prev_tsec,&t1,&t2](float p)
            {
                clock_gettime(CLOCK_MONOTONIC,&t2);
                size_t nsecs = 1000000000uLL*(t2.tv_sec-t1.tv_sec)
                    +(t2.tv_nsec-t1.tv_nsec);
                size_t secs = nsecs/1000000000;
                i32 percent = (i32)(100.0*p);
                size_t tsec = t2.tv_sec;
                // output when % changes or every second
                if (percent > prev_percent || tsec > prev_tsec)
                {
                    std::cerr << "\r" << "rendering... " << percent << "% ("
                        << secs << " sec elapsed)";
                    prev_percent = percent;
                    prev_tsec = tsec;
                }
            },
            [](const std::thread& thread, size_t index)
            {
                std::cerr << "starting thread " << index << " ("
                    << thread.get_id() << ")" << std::endl;
            });
        std::cerr << "\r" << "rendering... 100%" << std::endl;
        clock_gettime(CLOCK_MONOTONIC,&t2);
        size_t nsecs = 1000000000uLL*(t2.tv_sec-t1.tv_sec)
            +(t2.tv_nsec-t1.tv_nsec);
        float secs = (float)nsecs/1000000000.0;
        std::cerr << "render done" << std::endl;
        fprintf(stderr,"time (seconds): %f\n",secs);
        fprintf(stderr,"samples iterated: %lu\n",renderer.getSamplesIterated());
        fprintf(stderr,"samples plotted: %lu\n",renderer.getSamplesPlotted());
        fprintf(stderr,"samples/sec: %lf\n",renderer.getSamplesIterated()/secs);
        double ratio = (double)renderer.getSamplesPlotted()
            /renderer.getSamplesIterated();
        fprintf(stderr,"plotted/iterated: %lf%%\n",100.0*ratio);
        std::cerr << "xform selection:";
        for (size_t i = 0; i < flame.getXForms().size(); ++i)
            std::cerr << " " << renderer.getXFormFrequency()[i];
        std::cerr << std::endl;
        fprintf(stderr,"x min: %le\n",renderer.getPointExtremes()[0].first);
        fprintf(stderr,"x max: %le\n",renderer.getPointExtremes()[0].second);
        fprintf(stderr,"y min: %le\n",renderer.getPointExtremes()[1].first);
        fprintf(stderr,"y max: %le\n",renderer.getPointExtremes()[1].second);
        fprintf(stderr,"bad values: %lu\n",renderer.getBadValueCount());
        std::cerr << "bad value xforms:";
        for (size_t i = 0; i < renderer.getBadValueCount(); ++i)
            std::cerr << " " << renderer.getBadValueXForms()[i];
        std::cerr << std::endl;
        std::cerr << "bad value points:";
        for (auto p : renderer.getBadValuePoints())
            fprintf(stderr," (%le,%le)",p.x(),p.y());
        std::cerr << std::endl;
        size_t buffer_sum;
        hist_t sample_min,sample_max;
        renderer.histogramStats(buffer_sum,sample_min,sample_max);
        fprintf(stderr,"sample min: %u\n",sample_min);
        fprintf(stderr,"sample max: %u\n",sample_max);
        fprintf(stderr,"buffer sum: %lu\n",buffer_sum);
        size_t missed_samples = renderer.getSamplesPlotted()
            - (buffer_sum - buffer_sum_initial);
        fprintf(stderr,"missed samples: %lu\n",missed_samples);
#if 0
        // show histogram frequency information
        hist_t freq[8*sizeof(hist_t)]; // needing 1,2,3,..,32(or 64) bits
        for (size_t i = 0; i < 8*sizeof(hist_t); ++i)
            freq[i] = 0;
        hist_t *histptr = renderer.getHistogram();
        for (size_t i = 0; i < renderer.getHistogramSize(); ++i)
        {
            size_t j = 1;
            hist_t bit = 1;
            hist_t h = histptr[i];
            while (j < 8*sizeof(hist_t) && h >= bit)
            {
                ++j;
                bit <<= 1;
            }
            ++freq[j-1];
        }
        // output bit frequency with cumulative sum
        size_t cumulative = 0;
        for (size_t i = 0; i < 8*sizeof(hist_t); ++i)
        {
            cumulative += freq[i];
            fprintf(stderr,"freq(%lu bits) = %u (%f %%)\n",i+1,freq[i],
                100.0*((float)(cumulative)/renderer.getHistogramSize()));
        }
#endif
    }
    else
        std::cerr << "skipping render (0 samples)" << std::endl;
    std::cerr << "writing " << arg_output << std::endl;
    if (arg_type == "buf")
    {
        if (arg_output == "-")
        {
            if (!renderer.writeHistogram(std::cout))
            {
                std::cerr << "ERROR: failed writing to stdout" << std::endl;
                return 1;
            }
        }
        else
        {
            std::ofstream ofs(arg_output,std::ios::out|std::ios::binary);
            if (!ofs)
            {
                std::cerr << "ERROR: failed opening output file" << std::endl;
                return 1;
            }
            if (!renderer.writeHistogram(ofs))
            {
                std::cerr << "ERROR: failed writing output file" << std::endl;
                return 1;
            }
            ofs.close();
        }
        return 0;
    }
    std::ofstream ofs;
    if (arg_output != "-")
        ofs.open(arg_output,std::ios::out|std::ios::binary);
    std::ostream& os = arg_output != "-" ? ofs : std::cout;
    u8 *img8 = nullptr;
    u16 *img16 = nullptr;
    size_t X = renderer.getFlame().getSize()[0];
    size_t Y = renderer.getFlame().getSize()[1];
    std::function<num_t(hist_t)> scale;
    if (arg_scaler == "binary")
        scale = [](hist_t n) { return n ? 1.0 : 0.0; };
    else if (arg_scaler == "linear")
        scale = [](hist_t n) { return (num_t)n; };
    else if (arg_scaler == "log")
        scale = [](hist_t n) { return log(1.0+(num_t)n); };
    bool success;
    if (arg_img_bits == 8)
    {
        img8 = renderer.renderImageBuffer<u8>(scale);
        if (arg_type == "pgm")
            success = write_pgm(os,X,Y,img8);
        else // png
            success = write_png(os,X,Y,img8);
    }
    else
    {
        img16 = renderer.renderImageBuffer<u16>(scale);
        if (arg_type == "pgm")
            success = write_pgm(os,X,Y,img16);
        else // png
            success = write_png(os,X,Y,img16);
    }
    if (!success)
    {
        std::cerr << "ERROR: cannot write image" << std::endl;
        return 1;
    }
    std::cerr << "output done" << std::endl;
    // clean up and exit
    if (arg_output != "-")
        ofs.close();
    if (img8)
        delete[] img8;
    if (img16)
        delete[] img16;
    return 0;
}
