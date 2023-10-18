#include <iostream>
#include <cstring>
#include "tgaimage.hpp"

TGAImage::TGAImage(const int width, const int height,const int bpp){
    this->__width=width;this->__height=height;
    this->bpp=bpp;
}

bool TGAImage::readTgaFile(const std::string filename){
    std::ifstream in;
    in.open(filename,std::ios::binary);
    if(!in.is_open()){
        std::cerr<<"can't open file :" <<filename<<std::endl;
        return false;
    }
    TGAHeader header;
    in.read(reinterpret_cast<char*>(&header),sizeof(TGAHeader));
    if (in.bad()) {
        std::cerr<<"read error occurred while reading the header"<<std::endl;
        return false;
    }
    __width=header.width;
    __height=header.height;
    bpp=header.bitsPerPixel;
    if (__width<=0||__height<=0||(bpp!=GRAYSCALE&&bpp!=RGB&&bpp!=RGBA)) {
        std::cerr<<"bad bpp(or width or height) value"<<std::endl;
        return false;
    }
    size_t bytes=bpp*__width*__height;
    data=std::vector<uint8_t>(bytes,0);
    if (header.dataTypeCode==3||header.dataTypeCode==2){
        in.read(reinterpret_cast<char*>(data.data()),bytes);
        if (in.bad()){
            std::cerr<<"an error occurred while reading data"<<std::endl;
            return false;
        }
    }else if (header.dataTypeCode==10||header.dataTypeCode==11){
        if (!load_rle_data(in)){
            std::cerr<<"an error occurred while reading data"<<std::endl;
        }
    }else {
        std::cerr<<"unkown file format"<<(int)header.dataTypeCode<<std::endl;
        return false;
    }
    if (!header.imageDescriptor&0x20){
        flip_vertically();
    }
    if (!header.imageDescriptor&0x10){
        flip_horizontally();
    }
    std::cerr<<__width<<"x"<<__height<<"/"<<bpp*8<<std::endl;
    return true;
}

bool TGAImage::load_rle_data(std::ifstream& in){
    size_t pixelCount=__width*__height;
    size_t curPixel=0;
    size_t currentByte=0;
    while (curPixel<pixelCount){
        uint8_t chunk_header=0;
        chunk_header=in.get();
        if (in.bad()){
            std::cerr<<"an error occurred while reading data"<<std::endl;
            return false;
        }
        if (chunk_header<128){
            chunk_header++;
            for (int i=0;i<chunk_header;i++){
                
            }
        }
    }
}

bool TGAImage::writeTgaFile(const std::string filename,const bool vflip,const bool rle) const{
    constexpr std::uint8_t developer_area_ref[]={0,0,0,0};
    constexpr std::uint8_t extension_area_ref[]={0,0,0,0};
    constexpr std::uint8_t footer[18]={'T','R','U','E','V','I','S','I','O','N','-','X','F','I','L','E','.','\0'};
    std::ofstream out;
    out.open(filename,std::ios::binary);
    if(!out.is_open()){
        std::cerr<<"can't open file :" <<filename<<std::endl;
        return false;
    }
    TGAHeader header;
    header.bitsPerPixel = bpp<<3;
    header.width = __width; header.height = __height; 
    header.dataTypeCode = (bpp==GRAYSCALE) ?0x00:0x20;
    out.write(reinterpret_cast<const char*>(&header),sizeof(header));
    if (!out.good()){
        std::cerr<<"write error occurred while writing the header"<<std::endl;
        return false;
    }
    if (!rle){
        out.write(reinterpret_cast<const char*>(data.data()),data.size()*sizeof(uint8_t));
        if (!out.good()){
            std::cerr<<"can't unload raw data"<<std::endl;
            return false;
        }
    }
    else if (!unload_rle_data(out)){
        std::cerr<<"can't unload rle data"<<std::endl;
        return false;
    }
    out.write(reinterpret_cast<const char*>(developer_area_ref),sizeof(developer_area_ref));
}