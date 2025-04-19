/*
Renderer for the image output
TODO this would be better as a GUI to adjust and see immedeiate results

options:
-h,--help         show usage help message
-f,--flame        flame json file (required)
-i,--input        input buffers to add (can use multiple)
-o,--output       output file (currently png only)
-y,--gamma        gamma value for brightness adjustment
-m,--monochrome   binary image from histogram
-g,--grayscale    grayscale image from histogram
-c,--color        color image (currently 3d rgb color only)
-b,--bits         bit depth for grayscale/color (8 or 16)

planned options:
-s,--scaler use alternatives besides log
-1,--palette 1d palette coloring
-2,--hsl 2d hsl coloring
-3,--rgb 3d rgb coloring
-4,--rgba 4d rgba coloring

planned ideas:
- (variable width) blur filter
- supersampling
- more coloring modes
- more scaling modes
*/

#include "renderers/image_renderer.hpp"

#include "utils/endian.hpp"
#include "utils/image.hpp"
#include "utils/json.hpp"

#include <boost/program_options.hpp>

#include <fstream>
#include <iostream>
#include <string>
#include <vector>

const std::string VERSION = "unspecified";

namespace bpo = boost::program_options;

// variables for command line arguments
std::string arg_flame;
std::string arg_output;
std::vector<std::string> arg_input;
tkoz::flame::num_t arg_gamma = 1.0;
size_t arg_bits = 8;
bool arg_m,arg_g,arg_c;

// parsed data
tkoz::flame::Json json_flame;

// forward declaration
bool render_image(std::ostream&,tkoz::flame::ImageRenderer<>&);

int main(int argc, char **argv)
{
    // setup the argument parsing with boost
    bpo::options_description options("ffr-img usage");
    options.add_options()
        ("help,h","show this help message")
        ("flame,f",bpo::value<std::string>()->required(),"flame parameters JSON file")
        ("input,i",bpo::value<std::vector<std::string>>(),"buffers to add")
        ("output,o",bpo::value<std::string>()->required(),"output file")
        ("gamma,y",bpo::value<tkoz::flame::num_t>()->default_value(arg_gamma),"gamma value")
        ("monochrome,m",bpo::bool_switch(),"binary image from histogram")
        ("grayscale,g",bpo::bool_switch(),"grayscale image from histogram")
        ("color,c",bpo::bool_switch(),"color image from 3d color only")
        ("bits,b",bpo::value<size_t>()->default_value(arg_bits),"bits per color channel");
    bpo::variables_map args_map;
    bpo::store(bpo::command_line_parser(argc,argv).options(options).run(),args_map);

    // show help message
    if (args_map.count("help") || args_map.empty() || argc < 2)
    {
        std::cerr << options;
        return 1;
    }

    // parse arguments
    arg_flame = args_map["flame"].as<std::string>();
    arg_output = args_map["output"].as<std::string>();
    if (args_map.find("input") != args_map.end())
        arg_input = args_map["input"].as<std::vector<std::string>>();
    arg_gamma = args_map["gamma"].as<tkoz::flame::num_t>();
    if (arg_gamma < tkoz::flame::eps_v<tkoz::flame::num_t>)
        throw std::runtime_error("gamma too small");
    arg_m = args_map["monochrome"].as<bool>();
    arg_g = args_map["grayscale"].as<bool>();
    arg_c = args_map["color"].as<bool>();
    arg_bits = args_map["bits"].as<size_t>();
    if (arg_bits != 8 && arg_bits != 16)
        throw std::runtime_error("bits per channel must be 8 or 16");

    // print options
    std::cerr << "ffr-img version " << VERSION << std::endl;
    std::cerr << "--flame " << arg_flame << std::endl;
    for (auto& s : arg_input)
        std::cerr << "--input " << s << std::endl;
    std::cerr << "--output " << arg_output << std::endl;
    std::cerr << "--gamma " << arg_gamma << std::endl;
    if (arg_m)
        std::cerr << "--monochrome" << std::endl;
    if (arg_g)
        std::cerr << "--grayscale" << std::endl;
    if (arg_c)
        std::cerr << "--color" << std::endl;
    std::cerr << "--bits " << arg_bits << std::endl;
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
    if (dims != 2)
    {
        std::cerr << "ERROR: only 2D flames supported" << std::endl;
        return 1;
    }

    // create image renderer
    tkoz::flame::Flame<2> flame(json_flame);
    tkoz::flame::ImageRenderer<> img_ren(flame);
    auto& buf_ren = img_ren.getBufferRenderer();
    std::cerr << "buffer: " << (buf_ren.getBuffer().size()
        *sizeof(tkoz::flame::hist_t)) << " bytes, ";
#if TKOZ_LITTLE_ENDIAN
    std::cerr << "little endian, ";
#else
    std::cerr << "big endian, ";
#endif
    std::cerr << sizeof(tkoz::flame::hist_t) << " byte numbers" << std::endl;
    size_t color_dims = buf_ren.getColorDims();
    std::cerr << "color: " << color_dims << " dimensions" << std::endl;

    // check color type
    size_t color_flag_count = 0;
    if (arg_m) ++color_flag_count;
    if (arg_g) ++color_flag_count;
    if (arg_c) ++color_flag_count;
    if (color_flag_count != 1)
    {
        std::cerr << "ERROR: must choose exactly 1 coloring flag "
            << "(--monochrome, --grayscale, --color)" << std::endl;
        return 1;
    }

    // add input buffers
    if (arg_input.empty())
    {
        std::cerr << "ERROR: no input buffer files specified" << std::endl;
        return 1;
    }
    for (auto& s : arg_input)
    {
        std::cerr << "adding input file " << s << std::endl;
        if (s == "-")
        {
            if (!buf_ren.addBuffer(std::cin))
                throw std::runtime_error("error reading file");
        }
        else
        {
            std::ifstream input_file(s);
            if (!buf_ren.addBuffer(input_file))
                throw std::runtime_error("error reading file");
            input_file.close();
        }
    }

    // render image and output
    std::cerr << "processing buffer" << std::endl;
    bool write_success;
    if (arg_output == "-")
        write_success = render_image(std::cout,img_ren);
    else
    {
        std::ofstream output_file(arg_output);
        write_success = render_image(output_file,img_ren);
        output_file.close();
    }
    if (!write_success)
        throw std::runtime_error("error writing output");
    return 0;
}

bool render_image(std::ostream& os, tkoz::flame::ImageRenderer<>& img_ren)
{
    using num_t = tkoz::flame::num_t;
    using hist_t = tkoz::flame::hist_t;
    using u8 = tkoz::flame::u8;
    using u16 = tkoz::flame::u16;
    using buf_elem_t = tkoz::flame::buf_elem_t;
    using num_t3 = std::tuple<num_t,num_t,num_t>;
    //using cf_bool = tkoz::flame::cell_func_t<bool>;
    using cf_hist_t = tkoz::flame::util::cell_func_t<hist_t>;
    using cf_num_t = tkoz::flame::util::cell_func_t<num_t>;
    using cf_num_t3 = tkoz::flame::util::cell_func_t<num_t3>;

    // histogram min/max calculation
    cf_hist_t funch = [](const buf_elem_t *cell, size_t r) -> hist_t
    {
        (void)r;
        return cell->uintval;
    };
    auto [minh,maxh] = img_ren.getValueBounds(funch);
    std::cerr << "histogram bounds: " << minh << " " << maxh << std::endl;

    // log scaler
    cf_num_t func1 = [](const buf_elem_t *cell, size_t r) -> num_t
    {
        (void)r;
        return std::log(1 + (num_t)(cell->uintval));
    };

    // compute max with log scaling
    auto [min,max] = img_ren.getValueBounds(func1);
    std::cerr << "scaler bounds: " << min << " " << max << std::endl;
    if (max < tkoz::flame::eps_v<num_t>)
        throw std::runtime_error("histogram is (probably) empty");

    // gamma value
    num_t gp = 1.0 / arg_gamma;

    // log scaler to range [0,1]
    cf_num_t func2 = [max,gp](const buf_elem_t *cell, size_t _r) -> num_t
    {
        (void)_r;
        num_t l = std::log(1 + (num_t)(cell->uintval)) / max;
        return std::pow(l,gp);
    };

    if (arg_m)
    {
        /*
        This doesn't work because boost won't compile with gray1_image_t
        Using gray8_image_t to handle bitmap output for now
        // true for nonzero counter
        cf_bool func = [](const buf_elem_t *cell, size_t _r) -> bool
        {
            (void)_r;
            return cell->uintval != 0;
        };
        auto img = img_ren.renderBinaryImage(func);
        tkoz::flame::writePng(img,os);
        */
        cf_num_t func = [](const buf_elem_t *cell, size_t r) -> num_t
        {
            (void)r;
            return cell->uintval != 0 ? 1.0 : 0.0;
        };
        auto img = img_ren.renderGrayImage<u8>(func);
        tkoz::flame::writePng(img,os);
    }
    else if (arg_g)
    {
        if (arg_bits == 8)
        {
            auto img = img_ren.renderGrayImage<u8>(func2);
            tkoz::flame::writePng(img,os);
        }
        else
        {
            auto img = img_ren.renderGrayImage<u16>(func2);
            tkoz::flame::writePng(img,os);
        }
    }
    else if (arg_c)
    {
        if (img_ren.getColorDims() != 3)
            throw std::runtime_error("buffer must use 3 color dimensions");
        cf_num_t3 func = [max,gp](const buf_elem_t *cell, size_t _r) -> num_t3
        {
            (void)_r;
            num_t h = cell->uintval;
            num_t l = std::log(1 + (num_t)(cell->uintval)) / max;
            num_t ll = std::pow(l,gp);
            num_t r = cell[1].floatval / h;
            num_t g = cell[2].floatval / h;
            num_t b = cell[3].floatval / h;
            return {ll*r,ll*g,ll*b};
        };
        if (arg_bits == 8)
        {
            auto img = img_ren.renderColorImageRGB<u8>(func);
            tkoz::flame::writePng(img,os);
        }
        else
        {
            auto img = img_ren.renderColorImageRGB<u16>(func);
            tkoz::flame::writePng(img,os);
        }
    }
    else
        throw std::runtime_error("no coloring flag");
    return os.good();
}
