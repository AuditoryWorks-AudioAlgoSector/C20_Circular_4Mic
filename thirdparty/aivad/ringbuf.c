#include <stdlib.h>
#include <string.h>
#include "ringbuf.h"

/**
 * ringbuffer 简化版本
 * 写入每次写入一帧的目标数据
 * 读取则是每次都读取全部内容
*/

void rb_init(rb_t *rb, int buff_size) {

    rb->buffer = malloc(sizeof(float) * buff_size);
    memset(rb->buffer, 0, sizeof(float) * buff_size);

    rb->read_offset = 0;
    rb->write_offset = 0;
    rb->valid_size = 0;
    rb->total_size = buff_size;
}


void rb_deinit(rb_t *rb) {
    if (rb->buffer != NULL) {
        free(rb->buffer);
    }
}


void rb_write(rb_t *rb, float *buffer_to_write, int size) {

    memcpy(rb->buffer + rb->write_offset, buffer_to_write, size * sizeof(float));

    rb->write_offset += size;
    rb->write_offset %= rb->total_size;

    if(rb->valid_size < rb->total_size){
        rb->valid_size += size;
    }else{
        rb->valid_size = rb->total_size;
    }

}


void rb_read(rb_t *rb, float *buff_to_read, int size) {

    int first_read_size = 0;
    int offset = rb->write_offset;
    first_read_size = rb->total_size - offset;
    if (rb->write_offset==0){
        memcpy(buff_to_read, rb->buffer , sizeof(float) * first_read_size);
    } else{
        memcpy(buff_to_read, rb->buffer + offset, sizeof(float )*first_read_size);
        memcpy(buff_to_read + first_read_size, rb->buffer, sizeof(float)*(size - first_read_size));
    }
}

void rb_s16_init(rb_s16_t *rb_s16, int buff_size) {

    rb_s16->buffer = malloc(sizeof(short) * buff_size);
    memset(rb_s16->buffer, 0, sizeof(short) * buff_size);

    rb_s16->read_offset = 0;
    rb_s16->write_offset = 0;
    rb_s16->valid_size = 0;
    rb_s16->total_size = buff_size;
}


void rb_s16_deinit(rb_s16_t *rb_s16) {
    if (rb_s16->buffer != NULL) {
        free(rb_s16->buffer);
    }
}


void rb_s16_write(rb_s16_t *rb_s16, short *buffer_to_write, int size) {

    memcpy(rb_s16->buffer + rb_s16->write_offset, buffer_to_write, size * sizeof(short));

    rb_s16->write_offset += size;
    rb_s16->write_offset %= rb_s16->total_size;

}


void rb_s16_read(rb_s16_t *rb_s16, short *buff_to_read, int size) {
    int offset = rb_s16->write_offset;
    memcpy(buff_to_read, rb_s16->buffer+offset, sizeof(short) * size);

}
