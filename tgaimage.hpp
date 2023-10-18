#pragma once
#include <cstdint>
#include <fstream>
#include <vector>

#pragma pack(push,1)

struct TGAHeader{
    
    std::uint8_t idLength=0;
    std::uint8_t colorMapType=0;
    std::uint8_t dataTypeCode=0;
    std::uint8_t colorMapDepth=0;
    std::uint16_t colorMapOrigin=0;
    std::uint16_t colorMapLength=0;
    std::uint16_t x_Origin=0;
    std::uint16_t y_Origin=0;
    std::uint16_t width=0;
    std::uint16_t height=0;
    std::uint8_t bitsPerPixel=0;
    std::uint8_t imageDescriptor=0;
};
#pragma pack(pop)

struct TGAColor {
    std::uint8_t bgra[4]={0,0,0,0};
    std::uint8_t byteSpp=4;
    std::uint8_t& operator[] (const int i) {
        return bgra[i];
    }
};

class TGAImage {
bool load_rle_data(std::ifstream& in);
bool unload_rle_data(std::ofstream& out) const;
int __width;
int __height;
std::uint8_t bpp=0;
std::vector<std::uint8_t> data={};
public:
    enum Format {GRAYSCALE=1,RGB=3,RGBA=4};
    TGAImage()=default;
    TGAImage(const int width, const int height,const int hpp);
    bool readTgaFile(const std::string filename);
    bool writeTgaFile(const std::string filename,const bool vFlip=true,const bool rle=true) const;
    void flip_horizontally();
    void flip_vertically();
    TGAColor get(const int x,const int y) const;
    void set(const int x,const int y,const TGAColor& color);
    int width() const;
    int height() const;
};