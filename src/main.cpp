/*
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

planned options (not available yet):
[-p --precision]: calculation precision (single or double) (default single)
[-w --hist_bits]: histogram integer bit size (32 or 64) (default 32)
[-R --rng]: random number generator (java,isaac32,isaac64) (default isaac32)
[-r --seed]: random number generator seed seed (default random)
*/

#include <ctime>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>

#include <boost/program_options.hpp>

#include "renderers/renderer.hpp"
#include "utils/image.hpp"
#include "utils/misc.hpp"

const std::string VERSION = "unspecified";

namespace bpo = boost::program_options;
typedef float num_t; // can also be double (slower)

int main(int argc, char **argv)
{
    bpo::options_description options("ffbuf usage");
    options.add_options()
        ("help,h","show this message")
        ("flame,f",bpo::value<std::string>(),
            "flame parameters JSON file (required)")
        ("output,o",bpo::value<std::string>(),
            "output file (required)")
        ("input,i",bpo::value<std::string>()->default_value(""),
            "buffer (default none, render new buffer)")
        ("samples,s",bpo::value<size_t>()->default_value(0),
            "samples to render (default 0)")
        ("type,t",bpo::value<std::string>()->default_value(""),
            "output type (png,pgm,buf) (default use file extension)")
        ("img_bits,b",bpo::value<size_t>()->default_value(8),
            "bit depth for png or pgm output (8 or 16) (default 8)")
        ("threads,T",bpo::value<size_t>()->default_value(1),
            "number of threads to use (default 1)")
        ("batch_size,z",bpo::value<size_t>()->default_value(65536),
            "multithreading batch size (default 65536)")
        ("bad_values,B",bpo::value<size_t>()->default_value(10),
            "bad value limit for terminating render (default 10)")
        ("scaler,m",bpo::value<std::string>()->default_value("log"),
            "scaling function for image render (bin,lin,log) (default log)");
    bpo::variables_map args;
    bpo::store(bpo::command_line_parser(argc,argv).options(options).run(),args);
    if (args.count("help") || args.empty())
    {
        std::cerr << options;
        return 1;
    }
    if (!args.count("flame")) // required argument
    {
        std::cerr << "error: flame JSON not specified" << std::endl;
        return 1;
    }
    if (!args.count("output")) // required argument
    {
        std::cerr << "error: output file not specified" << std::endl;
        return 1;
    }
    // cmdline values
    std::string arg_flame = args["flame"].as<std::string>();
    std::string arg_input = args["input"].as<std::string>();
    size_t arg_samples = args["samples"].as<size_t>();
    std::string arg_output = args["output"].as<std::string>();
    std::string arg_type = args["type"].as<std::string>();
    size_t arg_img_bits = args["img_bits"].as<size_t>();
    size_t arg_threads = args["threads"].as<size_t>();
    size_t arg_batch_size = args["batch_size"].as<size_t>();
    size_t arg_bad_values = args["bad_values"].as<size_t>();
    std::string arg_scaler = args["scaler"].as<std::string>();
    //std::string arg_precision = args["precision"].as<std::string>();
    //size_t arg_hist_bits = args["hist_bits"].as<size_t>();
    // check args
    if (arg_type != "" && arg_type != "png" && arg_type != "pgm"
        && arg_type != "buf")
    {
        std::cerr << "error: invalid output type" << std::endl;
        return 1;
    }
    if (arg_img_bits != 8 && arg_img_bits != 16)
    {
        std::cerr << "error: invalid image bits per pixel" << std::endl;
        return 1;
    }
    if (arg_threads < 1 || arg_threads > 128)
    {
        std::cerr << "error: must use between 1 and 128 threads" << std::endl;
        return 1;
    }
    if (arg_batch_size < (1<<8) || arg_batch_size > (1<<30))
    {
        std::cerr << "error: batch size < 2^8 or > 2^30" << std::endl;
        return 1;
    }
    if (arg_bad_values < 1)
    {
        std::cerr << "error: bad values must be positive" << std::endl;
        return 1;
    }
    if (arg_scaler != "binary" && arg_scaler != "linear" && arg_scaler != "log")
    {
        std::cerr << "error: scaler must be binary/linear/log" << std::endl;
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
            std::cerr << "error: unable to infer output type" << std::endl;
            return 1;
        }
    }
    std::cerr << "=== FFBUF version " << VERSION << " ===" << std::endl;
    // print options
    std::cerr << "command line arguments:" << std::endl;
    std::cerr << "ffbuf" << std::endl;
    std::cerr << "--flame: " << arg_flame << std::endl;
    std::cerr << "--input: " << arg_input << std::endl;
    std::cerr << "--samples: " << arg_samples << std::endl;
    std::cerr << "--output: " << arg_output << std::endl;
    std::cerr << "--type: " << arg_type << std::endl;
    std::cerr << "--img_bits: " << arg_img_bits << std::endl;
    std::cerr << "--threads: " << arg_threads << std::endl;
    std::cerr << "--batch_size: " << arg_batch_size << std::endl;
    std::cerr << "--bad_values: " << arg_bad_values << std::endl;
    std::cerr << "--" << std::endl;
    // parse flame file
    Json json_flame;
    if (arg_flame == "-") // input from stdin
        json_flame = Json(std::cin);
    else
        json_flame = Json(read_text_file(arg_flame));
    std::cerr << "flame json (comments removed): " << json_flame << std::endl;
    tkoz::flame::RendererBasic<num_t,u32> renderer(json_flame);
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
    u32 *buf = renderer.getHistogram(); // buffer to overwrite
    // load buffer if specified
    if (arg_input != "")
    {
        if (arg_input == "-") // input from stdin
        {
            if (!std::cin.read((char*)buf,renderer.getHistogramSizeBytes()))
            {
                std::cerr << "error: unable to read (enough) buffer bytes"
                    << std::endl;
                return 1;
            }
            if (!std::cin.eof())
                std::cerr << "warn: did not reach EOF of stdin" << std::endl;
        }
        else
        {
            std::ifstream ifs(arg_input,std::ios::in|std::ios::binary);
            if (!ifs.read((char*)buf,renderer.getHistogramSizeBytes()))
            {
                std::cerr << "error: unable to read (enough) buffer bytes"
                    << std::endl;
                return 1;
            }
            if (!ifs.eof())
                std::cerr << "warn: did not reach EOF of file" << std::endl;
            ifs.close();
        }
        fprintf(stderr,"read %lu bytes\n",renderer.getHistogramSizeBytes());
    }
    size_t buffer_sum_initial = 0;
    for (size_t i = 0; i < renderer.getHistogramSize(); ++i)
        buffer_sum_initial += buf[i];
    fprintf(stderr,"buffer sum (initial): %lu\n",buffer_sum_initial);
    if (arg_samples)
    {
        std::cerr << "render start" << std::endl;
        struct timespec t1,t2;
        clock_gettime(CLOCK_MONOTONIC,&t1);
        i32 prev_percent = -1;
        size_t prev_tsec = t1.tv_sec;
        //renderer.renderBuffer(arg_samples,rng);
        renderer.renderBufferParallel(arg_samples,arg_threads,
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
            [](std::thread& thread, size_t index)
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
        for (size_t i = 0; i < renderer.getXFormsLength(); ++i)
            std::cerr << " " << renderer.getXFormDistribution()[i];
        std::cerr << std::endl;
        fprintf(stderr,"x min: %le\n",renderer.getXMin());
        fprintf(stderr,"x max: %le\n",renderer.getXMax());
        fprintf(stderr,"y min: %le\n",renderer.getYMin());
        fprintf(stderr,"y max: %le\n",renderer.getYMax());
        fprintf(stderr,"bad values: %lu\n",renderer.getBadValueCount());
        std::cerr << "bad value xforms:";
        for (size_t i = 0; i < renderer.getBadValueCount(); ++i)
            std::cerr << " " << renderer.getBadValueXForms()[i];
        std::cerr << std::endl;
        std::cerr << "bad value points:";
        for (auto p : renderer.getBadValuePoints())
            fprintf(stderr," (%le,%le)",p.x(),p.y());
        std::cerr << std::endl;
        u32 sample_min = -1;
        u32 sample_max = 0;
        size_t buffer_sum = 0;
        for (size_t i = 0; i < renderer.getHistogramSize(); ++i)
        {
            sample_min = std::min(sample_min,buf[i]);
            sample_max = std::max(sample_max,buf[i]);
            buffer_sum += buf[i];
        }
        fprintf(stderr,"sample min: %u\n",sample_min);
        fprintf(stderr,"sample max: %u\n",sample_max);
        fprintf(stderr,"buffer sum: %lu\n",buffer_sum);
        size_t missed_samples = renderer.getSamplesPlotted()
            - (buffer_sum - buffer_sum_initial);
        fprintf(stderr,"missed samples: %lu\n",missed_samples);
    }
    else
        std::cerr << "skipping render (0 samples)" << std::endl;
    std::cerr << "writing " << arg_output << std::endl;
    if (arg_type == "buf")
    {
        if (arg_output == "-")
        {
            std::cout.write((char*)buf,renderer.getHistogramSizeBytes());
            if (!std::cout)
            {
                std::cerr << "error: cannot write to stdout" << std::endl;
                return 1;
            }
        }
        else
        {
            std::ofstream ofs(arg_output,std::ios::out|std::ios::binary);
            ofs.write((char*)buf,renderer.getHistogramSizeBytes());
            if (!ofs)
            {
                std::cerr << "error: cannot write output file" << std::endl;
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
    std::function<num_t(u32)> scale;
    if (arg_scaler == "binary")
        scale = [](u32 n) { return n ? 1.0 : 0.0; };
    else if (arg_scaler == "linear")
        scale = [](u32 n) { return (num_t)n; };
    else if (arg_scaler == "log")
        scale = [](u32 n) { return log(1.0+(num_t)n); };
    bool success;
    if (arg_img_bits == 8)
    {
        img8 = renderer.renderImage<u8>(scale);
        if (arg_type == "pgm")
            success = write_pgm(os,X,Y,img8);
        else // png
            success = write_png(os,X,Y,img8);
    }
    else
    {
        img16 = renderer.renderImage<u16>(scale);
        if (arg_type == "pgm")
            success = write_pgm(os,X,Y,img16);
        else // png
            success = write_png(os,X,Y,img16);
    }
    if (!success)
    {
        std::cerr << "error: cannot write image" << std::endl;
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
