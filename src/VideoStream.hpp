#pragma once

#include <GL/gl.h>
#include <vector>

class VideoStream {
public:
    VideoStream(const char* filename);
    ~VideoStream();

    void writeFrame(const std::vector<uint8_t>& pixels);

private:
    struct VideoStream_impl *impl;
};
