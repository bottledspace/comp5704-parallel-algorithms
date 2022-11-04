#include "VideoStream.hpp"
#include <cassert>
#include <algorithm>

static void add_stream(OutputStream *ost, AVFormatContext *oc,
                       const AVCodec **codec,
                       enum AVCodecID codec_id)
{
    const AVChannelLayout layout = AV_CHANNEL_LAYOUT_STEREO;
    AVCodecContext *c;

    assert(*codec = avcodec_find_encoder(codec_id));
    assert(ost->tmp_pkt = av_packet_alloc());
    assert(ost->st = avformat_new_stream(oc, NULL));
    ost->st->id = oc->nb_streams-1;
    assert(c = avcodec_alloc_context3(*codec));
    ost->enc = c;

    c->codec_id = codec_id;
    c->bit_rate = 400000;
    c->width    = 512; // Resolution must be a multiple of two.
    c->height   = 512;
    /* timebase: This is the fundamental unit of time (in seconds) in terms
        * of which frame timestamps are represented. For fixed-fps content,
        * timebase should be 1/framerate and timestamp increments should be
        * identical to 1. */
    ost->st->time_base = (AVRational){ 1, STREAM_FRAME_RATE };
    c->time_base       = ost->st->time_base;
    c->gop_size      = 12; /* emit one intra frame every twelve frames at most */
    c->pix_fmt       = STREAM_PIX_FMT;
    if (c->codec_id == AV_CODEC_ID_MPEG2VIDEO) {
        /* just for testing, we also add B-frames */
        c->max_b_frames = 2;
    }
    if (c->codec_id == AV_CODEC_ID_MPEG1VIDEO) {
        /* Needed to avoid using macroblocks in which some coeffs overflow.
            * This does not happen with normal video, it just happens here as
            * the motion of the chroma plane does not match the luma plane. */
        c->mb_decision = 2;
    }
    /* Some formats want stream headers to be separate. */
    if (oc->oformat->flags & AVFMT_GLOBALHEADER)
        c->flags |= AV_CODEC_FLAG_GLOBAL_HEADER;
}

static AVFrame *alloc_picture(enum AVPixelFormat pix_fmt, int width, int height)
{
    AVFrame *picture;
    assert(picture = av_frame_alloc());
    picture->format = pix_fmt;
    picture->width  = width;
    picture->height = height;
    assert(av_frame_get_buffer(picture, 0) >= 0);
    return picture;
}

static void open_video(AVFormatContext *oc, const AVCodec *codec,
                       OutputStream *ost, AVDictionary *opt_arg)
{
    AVCodecContext *c = ost->enc;
    AVDictionary *opt = NULL;

    av_dict_copy(&opt, opt_arg, 0);
    assert(avcodec_open2(c, codec, &opt) >= 0);
    av_dict_free(&opt);
    ost->frame = alloc_picture(c->pix_fmt, c->width, c->height);
    assert(ost->frame);
    ost->tmp_frame = alloc_picture(AV_PIX_FMT_RGB32, c->width, c->height);
    assert(ost->tmp_frame);
    assert(avcodec_parameters_from_context(ost->st->codecpar, c) >= 0);
}

static void fill_yuv_image(AVFrame *pict, int frame_index, int width, int height)
{
    int x, y;
    glPixelStorei(GL_PACK_ROW_LENGTH, pict->linesize[0]/4);
    glReadPixels(0,0, width,height, GL_RGBA, GL_UNSIGNED_BYTE, pict->data[0]);
}

AVFrame *VideoStream::sendFrame(void)
{
    OutputStream *ost = &video_st;
    AVCodecContext *c = ost->enc;

    /* check if we want to generate more frames */
    /*if (av_compare_ts(ost->next_pts, c->time_base,
                      STREAM_DURATION, (AVRational){ 1, 1 }) > 0)
        return NULL;*/
    assert(av_frame_make_writable(ost->frame) >= 0);

    if (!ost->sws_ctx) {
        ost->sws_ctx = sws_getContext(c->width, c->height,
                                        AV_PIX_FMT_RGB32,
                                        c->width, c->height,
                                        c->pix_fmt,
                                        SCALE_FLAGS, NULL, NULL, NULL);
        assert(ost->sws_ctx);
    }
    fill_yuv_image(ost->tmp_frame, ost->next_pts, c->width, c->height);
    uint8_t* p = ost->tmp_frame->data[0];
    long rl = ost->tmp_frame->linesize[0];
    for (int y = 0; y < c->height/2; y++)
        std::swap_ranges(p+rl*y,p+rl*(y+1),p+rl*(c->height-y-1));
    sws_scale(ost->sws_ctx, (const uint8_t * const *)ost->tmp_frame->data,
                ost->tmp_frame->linesize, 0, c->height, ost->frame->data,
                ost->frame->linesize);

    ost->frame->pts = ost->next_pts++;
    return ost->frame;
}


bool VideoStream::writeFrame(void)
{
    int ret;
    AVFrame* frame = sendFrame();
    assert(avcodec_send_frame(video_st.enc, frame) >= 0);
    while (ret >= 0) {
        ret = avcodec_receive_packet(video_st.enc, video_st.tmp_pkt);
        if (ret == AVERROR(EAGAIN) || ret == AVERROR_EOF)
            break;
        else if (ret < 0) {
            fprintf(stderr, "Error encoding a frame: %s\n", "??");
            exit(1);
        }
        av_packet_rescale_ts(video_st.tmp_pkt, video_st.enc->time_base, video_st.st->time_base);
        video_st.tmp_pkt->stream_index = video_st.st->index;
        assert(av_interleaved_write_frame(oc, video_st.tmp_pkt) >= 0);
    }
    return ret == AVERROR_EOF ? 1 : 0;
}

VideoStream::~VideoStream()
{
    av_write_trailer(oc);
    avcodec_free_context(&video_st.enc);
    av_frame_free(&video_st.frame);
    av_frame_free(&video_st.tmp_frame);
    av_packet_free(&video_st.tmp_pkt);
    sws_freeContext(video_st.sws_ctx);
    swr_free(&video_st.swr_ctx);
    if (!(fmt->flags & AVFMT_NOFILE))
        avio_closep(&oc->pb);
    avformat_free_context(oc);
}

VideoStream::VideoStream(const char* filename)
{
    AVDictionary *opt = NULL;
    avformat_alloc_output_context2(&oc, NULL, NULL, filename);
    assert(oc);
    fmt = oc->oformat;
    add_stream(&video_st, oc, &video_codec, fmt->video_codec);
    open_video(oc, video_codec, &video_st, opt);
    av_dump_format(oc, 0, filename, 1);

    if (!(fmt->flags & AVFMT_NOFILE))
        assert(avio_open(&oc->pb, filename, AVIO_FLAG_WRITE) >= 0);
    assert(avformat_write_header(oc, &opt) >= 0);
}
