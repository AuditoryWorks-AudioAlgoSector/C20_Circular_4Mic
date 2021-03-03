#ifndef RINGBUF_H
#define RINGBUF_H
typedef struct rb {
    float *buffer;      //存放实际的数据
    int read_offset;  //读取地址相对buffer的偏移量
    int write_offset; //写入地址相对buffer的偏移量
    int valid_size;   //buffer的有效size
    int total_size;   //buffer的总大小
} rb_t;
typedef struct rb_s16 {
    short *buffer;      //存放实际的数据
    int read_offset;  //读取地址相对buffer的偏移量
    int write_offset; //写入地址相对buffer的偏移量
    int valid_size;   //buffer的有效size
    int total_size;   //buffer的总大小
} rb_s16_t;


#ifdef __cplusplus
extern "C"{
#endif
void rb_init(rb_t *rb, int buff_size);

void rb_deinit(rb_t *rb);

void rb_write(rb_t *rb, float *buffer_to_write, int size );

void rb_read(rb_t *rb, float *buff_to_read, int size);

void rb_s16_init(rb_s16_t *rb_s16, int buff_size);
void rb_s16_deinit(rb_s16_t *rb_s16);
void rb_s16_write(rb_s16_t *rb_s16, short *buffer_to_write, int size);
void rb_s16_read(rb_s16_t *rb_s16, short *buff_to_read, int size);

#ifdef __cplusplus
};
#endif

#endif//RINGBUF_H