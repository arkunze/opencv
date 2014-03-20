/*M///////////////////////////////////////////////////////////////////////////////////////
//
//  IMPORTANT: READ BEFORE DOWNLOADING, COPYING, INSTALLING OR USING.
//
//  By downloading, copying, installing or using the software you agree to this license.
//  If you do not agree to this license, do not download, install,
//  copy or use the software.
//
//
//                           License Agreement
//                For Open Source Computer Vision Library
//
// Copyright (C) 2010-2013, Advanced Micro Devices, Inc., all rights reserved.
// Copyright (C) 2014, Intel Corporation, all rights reserved.
// Third party copyrights are property of their respective owners.
//
// Redistribution and use in source and binary forms, with or without modification,
// are permitted provided that the following conditions are met:
//
//   * Redistribution's of source code must retain the above copyright notice,
//     this list of conditions and the following disclaimer.
//
//   * Redistribution's in binary form must reproduce the above copyright notice,
//     this list of conditions and the following disclaimer in the documentation
//     and/or other materials provided with the distribution.
//
//   * The name of the copyright holders may not be used to endorse or promote products
//     derived from this software without specific prior written permission.
//
// This software is provided by the copyright holders and contributors as is and
// any express or implied warranties, including, but not limited to, the implied
// warranties of merchantability and fitness for a particular purpose are disclaimed.
// In no event shall the Intel Corporation or contributors be liable for any direct,
// indirect, incidental, special, exemplary, or consequential damages
// (including, but not limited to, procurement of substitute goods or services;
// loss of use, data, or profits; or business interruption) however caused
// and on any theory of liability, whether in contract, strict liability,
// or tort (including negligence or otherwise) arising in any way out of
// the use of this software, even if advised of the possibility of such damage.
//
//M*/

#ifdef BORDER_REPLICATE
//BORDER_REPLICATE:     aaaaaa|abcdefgh|hhhhhhh
#define ADDR_L(i, l_edge, r_edge)  ((i) <  (l_edge) ? (l_edge)   : (i))
#define ADDR_R(i, r_edge, addr)    ((i) >= (r_edge) ? (r_edge)-1 : (addr))
#define ADDR_H(i, t_edge, b_edge)  ((i) <  (t_edge) ? (t_edge)   :(i))
#define ADDR_B(i, b_edge, addr)    ((i) >= (b_edge) ? (b_edge)-1 :(addr))
#endif

#ifdef BORDER_REFLECT
//BORDER_REFLECT:       fedcba|abcdefgh|hgfedcb
#define ADDR_L(i, l_edge, r_edge)  ((i) <  (l_edge) ? -(i)-1               : (i))
#define ADDR_R(i, r_edge, addr)    ((i) >= (r_edge) ? -(i)-1+((r_edge)<<1) : (addr))
#define ADDR_H(i, t_edge, b_edge)  ((i) <  (t_edge) ? -(i)-1 : (i))
#define ADDR_B(i, b_edge, addr)    ((i) >= (b_edge) ? -(i)-1+((b_edge)<<1) : (addr))
#endif

#ifdef BORDER_REFLECT_101
//BORDER_REFLECT_101:   gfedcb|abcdefgh|gfedcba
#define ADDR_L(i, l_edge, r_edge)  ((i) <  (l_edge) ? -(i)                 : (i))
#define ADDR_R(i, r_edge, addr)    ((i) >= (r_edge) ? -(i)-2+((r_edge)<<1) : (addr))
#define ADDR_H(i, t_edge, b_edge)  ((i) <  (t_edge) ? -(i)                 : (i))
#define ADDR_B(i, b_edge, addr)    ((i) >= (b_edge) ? -(i)-2+((b_edge)<<1) : (addr))
#endif

//blur function does not support BORDER_WRAP
#ifdef BORDER_WRAP
//BORDER_WRAP:          cdefgh|abcdefgh|abcdefg
#define ADDR_L(i, l_edge, r_edge)  ((i) <  (l_edge) ? (i)+(r_edge) : (i))
#define ADDR_R(i, r_edge, addr)    ((i) >= (r_edge) ? (i)-(r_edge) : (addr))
#define ADDR_H(i, t_edge, b_edge)  ((i) <  (t_edge) ? (i)+(b_edge) : (i))
#define ADDR_B(i, b_edge, addr)    ((i) >= (b_edge) ? (i)-(b_edge) : (addr))
#endif

#ifdef BORDER_ISOLATED
#define ISOLATED_MIN(VAL) (VAL)
#else
#define ISOLATED_MIN(VAL) 0
#endif

#define EXTRAPOLATE(x, y, minX, minY, maxX, maxY) \
    { \
        int _row = y - ISOLATED_MIN(minY), _col = x - ISOLATED_MIN(minX); \
        _row = ADDR_H(_row, 0, maxY - ISOLATED_MIN(minY)); \
        _row = ADDR_B(_row, maxY - ISOLATED_MIN(minY), _row); \
        y = _row + ISOLATED_MIN(minY); \
        \
        _col = ADDR_L(_col, 0, maxX - ISOLATED_MIN(minX)); \
        _col = ADDR_R(_col, maxX - ISOLATED_MIN(minX), _col); \
        x = _col + ISOLATED_MIN(minX); \
    }

#if USE_DOUBLE
#ifdef cl_amd_fp64
#pragma OPENCL EXTENSION cl_amd_fp64:enable
#elif defined (cl_khr_fp64)
#pragma OPENCL EXTENSION cl_khr_fp64:enable
#endif
#define FPTYPE double
#define CONVERT_TO_FPTYPE CAT(convert_double, VEC_SIZE)
#else
#define FPTYPE float
#define CONVERT_TO_FPTYPE CAT(convert_float, VEC_SIZE)
#endif

#if DATA_DEPTH == 0
#define BASE_TYPE uchar
#elif DATA_DEPTH == 1
#define BASE_TYPE char
#elif DATA_DEPTH == 2
#define BASE_TYPE ushort
#elif DATA_DEPTH == 3
#define BASE_TYPE short
#elif DATA_DEPTH == 4
#define BASE_TYPE int
#elif DATA_DEPTH == 5
#define BASE_TYPE float
#elif DATA_DEPTH == 6
#define BASE_TYPE double
#else
#error data_depth
#endif

#define __CAT(x, y) x##y
#define CAT(x, y) __CAT(x, y)

#define uchar1 uchar
#define char1 char
#define ushort1 ushort
#define short1 short
#define int1 int
#define float1 float
#define double1 double

#define convert_uchar1_sat_rte convert_uchar_sat_rte
#define convert_char1_sat_rte convert_char_sat_rte
#define convert_ushort1_sat_rte convert_ushort_sat_rte
#define convert_short1_sat_rte convert_short_sat_rte
#define convert_int1_sat_rte convert_int_sat_rte
#define convert_float1
#define convert_double1

#if DATA_DEPTH == 5 || DATA_DEPTH == 6
#define CONVERT_TO_TYPE CAT(CAT(convert_, BASE_TYPE), VEC_SIZE)
#else
#define CONVERT_TO_TYPE CAT(CAT(CAT(convert_, BASE_TYPE), VEC_SIZE), _sat_rte)
#endif

#define VEC_SIZE DATA_CHAN

#define VEC_TYPE CAT(BASE_TYPE, VEC_SIZE)
#define TYPE VEC_TYPE

#define SCALAR_TYPE CAT(FPTYPE, VEC_SIZE)

#define INTERMEDIATE_TYPE CAT(FPTYPE, VEC_SIZE)

#ifdef BORDER_CONSTANT
#define CONSTANT_BORDER_PARAM   SCALAR_TYPE borderValue,
#define CONSTANT_BORDER_ARG     borderValue
#else
#define CONSTANT_BORDER_PARAM
#define CONSTANT_BORDER_ARG     (SCALAR_TYPE)0
#endif

struct RectCoords
{
    int x1, y1, x2, y2;
};

//#define DEBUG
#ifdef DEBUG
#define DEBUG_ONLY(x) x
#define ASSERT(condition) do { if (!(condition)) { printf("BUG in boxFilter kernel (global=%d,%d): " #condition "\n", get_global_id(0), get_global_id(1)); } } while (0)
#else
#define DEBUG_ONLY(x) (void)0
#define ASSERT(condition) (void)0
#endif

#define vload1(OFFSET, PTR) (*(PTR + OFFSET))
#define PX_LOAD_VEC_TYPE CAT(BASE_TYPE, PX_LOAD_VEC_SIZE)
#define PX_LOAD_FLOAT_VEC_TYPE CAT(FPTYPE, PX_LOAD_VEC_SIZE)
#define PX_LOAD_FLOAT_VEC_CONV CAT(convert_, PX_LOAD_FLOAT_VEC_TYPE)
#define PX_LOAD CAT(vload, PX_LOAD_VEC_SIZE)

#ifdef BORDER_ISOLATED
inline bool isBorder(const struct RectCoords bounds, int2 coord, int numPixels)
{
    return (coord.x < bounds.x1 || coord.y < bounds.y1 || coord.x + numPixels > bounds.x2 || coord.y >= bounds.y2);
}
#else
inline bool isBorder(const struct RectCoords bounds, int2 coord, int numPixels)
{
    return (coord.x < 0 || coord.y < 0 || coord.x + numPixels > bounds.x2 || coord.y >= bounds.y2);
}
#endif

#ifdef BORDER_CONSTANT
INTERMEDIATE_TYPE getBorderPixel(const struct RectCoords bounds, int2 coord, SCALAR_TYPE constBorder,
                                 __global const uchar* srcptr, int srcstep)
{
    return CONVERT_TO_FPTYPE(constBorder);
}
#else
INTERMEDIATE_TYPE getBorderPixel(const struct RectCoords bounds, int2 coord, SCALAR_TYPE constBorder,
                                 __global const uchar* srcptr, int srcstep)
{
    int selected_col = coord.x;
    int selected_row = coord.y;

    EXTRAPOLATE(selected_col, selected_row,
            bounds.x1, bounds.y1,
            bounds.x2, bounds.y2
        );

    // debug border mapping
    //printf("pos=%d,%d --> %d, %d\n", pos.x, pos.y, selected_col, selected_row);

    coord = (int2)(selected_col, selected_row);
    if(coord.x >= 0 && coord.y >= 0 && coord.x < bounds.x2 && coord.y < bounds.y2)
    {
        __global VEC_TYPE* ptr = (__global VEC_TYPE*)(srcptr + mul24(coord.y, srcstep) + coord.x * sizeof(VEC_TYPE));
        return CONVERT_TO_FPTYPE(*ptr);
    }
    else
    {
        // for debug only
        DEBUG_ONLY(printf("BUG in filter2D kernel\n"));
        return (INTERMEDIATE_TYPE)(0.0f);
    }
}
#endif

inline INTERMEDIATE_TYPE readSrcPixelSingle(int2 pos, __global const uchar* srcptr,
                                            int srcstep, const struct RectCoords srcCoords,
                                            SCALAR_TYPE borderValue)
{
    if (!isBorder(srcCoords, pos, 1))
    {
        __global VEC_TYPE* ptr = (__global VEC_TYPE*)(srcptr + mul24(pos.y, srcstep) + pos.x * sizeof(VEC_TYPE));
        return CONVERT_TO_FPTYPE(*ptr);
    }
    else
    {
        //if (pos.x == 83 && pos.y == 4) printf("Input pixel @(%d, %d) is border\n", pos.x, pos.y);
        return getBorderPixel(srcCoords, pos, borderValue, srcptr, srcstep);
    }
}

inline PX_LOAD_FLOAT_VEC_TYPE readSrcPixelGroup(int2 pos, __global const uchar* srcptr,
                                                int srcstep, const struct RectCoords srcCoords)
{
    __global BASE_TYPE* ptr = (__global BASE_TYPE*)(srcptr + mul24(pos.y, srcstep) + pos.x * sizeof(VEC_TYPE));
    return PX_LOAD_FLOAT_VEC_CONV(PX_LOAD(0, ptr));
}

#define KERNEL_LOAD		CAT(vload, KERN_LOAD_VEC_SIZE)
#define KERNEL_LOAD_T	CAT(FPTYPE, KERN_LOAD_VEC_SIZE)

// Macros to ensure unrolled loops
#define LOOP1(VAR, STMT) (STMT); (VAR)++;
#define LOOP2(VAR, STMT) LOOP1(VAR, STMT); (STMT); (VAR)++;
#define LOOP3(VAR, STMT) LOOP2(VAR, STMT); (STMT); (VAR)++;
#define LOOP4(VAR, STMT) LOOP3(VAR, STMT); (STMT); (VAR)++;
#define LOOP5(VAR, STMT) LOOP4(VAR, STMT); (STMT); (VAR)++;
#define LOOP6(VAR, STMT) LOOP5(VAR, STMT); (STMT); (VAR)++;
#define LOOP7(VAR, STMT) LOOP6(VAR, STMT); (STMT); (VAR)++;
#define LOOP8(VAR, STMT) LOOP7(VAR, STMT); (STMT); (VAR)++;
#define LOOP9(VAR, STMT) LOOP8(VAR, STMT); (STMT); (VAR)++;
#define LOOP10(VAR, STMT) LOOP9(VAR, STMT); (STMT); (VAR)++;
#define LOOP11(VAR, STMT) LOOP10(VAR, STMT); (STMT); (VAR)++;
#define LOOP12(VAR, STMT) LOOP11(VAR, STMT); (STMT); (VAR)++;
#define LOOP13(VAR, STMT) LOOP12(VAR, STMT); (STMT); (VAR)++;

#define LOOP(N, VAR, STMT) CAT(LOOP, N)((VAR), (STMT))

__kernel
void filter2DSmall(__global const uchar* srcptr, int srcstep, int srcOffsetX, int srcOffsetY, int srcEndX, int srcEndY,
                __global uchar* dstptr, int dststep, int dstoffset,
               int rows, int cols,
               CONSTANT_BORDER_PARAM
               __constant FPTYPE* kernelData // transposed: [KERNEL_SIZE_X][KERNEL_SIZE_Y2_ALIGNED]
               )
{
    const struct RectCoords srcCoords = {srcOffsetX, srcOffsetY, srcEndX, srcEndY};

    const int startX = get_global_id(0) * PX_PER_WI_X;
    const int startY = get_global_id(1) * PX_PER_WI_Y;

    if ((startX >= cols) || (startY >= rows))
    {
        return;
    }

    FPTYPE privateKernel[KERNEL_SIZE_Y][KERNEL_SIZE_X];
    for (int kernY = 0; kernY < KERNEL_SIZE_Y; ++kernY)
    {
        for (int kernX = 0; kernX < KERNEL_SIZE_X; kernX += KERN_LOAD_VEC_SIZE)
        {
            KERNEL_LOAD_T kernData = KERNEL_LOAD(0, kernelData +
                                                 (kernY * KERNEL_SIZE_X) +
                                                 kernX);
            *(KERNEL_LOAD_T*)&privateKernel[kernY][kernX] = kernData;
        }
    }

    INTERMEDIATE_TYPE privateData[PX_PER_WI_Y + KERNEL_SIZE_Y - 1][PRIV_DATA_WIDTH];

    // Load all of the pixels needed for the calculation
    int py = 0;
    LOOP(PX_LOAD_Y_ITERATIONS, py,
    {
        int y = startY + py;
        int px = 0;
        LOOP(PX_LOAD_X_ITERATIONS, px,
        {
            int x = startX + (px * PX_LOAD_NUM_PX);
            int2 srcPos = (int2)(srcCoords.x1 + x - ANCHOR_X, srcCoords.y1 + y - ANCHOR_Y);

            if (!isBorder(srcCoords, srcPos, PX_LOAD_NUM_PX))
            {
                PX_LOAD_FLOAT_VEC_TYPE p = readSrcPixelGroup(srcPos, srcptr, srcstep, srcCoords);
                *((PX_LOAD_FLOAT_VEC_TYPE*)&privateData[py][px * PX_LOAD_NUM_PX]) = p;
            }
            else
            {
                int lx = 0;
                LOOP(PX_LOAD_NUM_PX, lx,
                {
                    INTERMEDIATE_TYPE p = readSrcPixelSingle(srcPos, srcptr, srcstep, srcCoords,
                                                             CONSTANT_BORDER_ARG);
                    *((INTERMEDIATE_TYPE*)&privateData[py][px * PX_LOAD_NUM_PX + lx]) = p;
                    srcPos.x++;
                });
            }
        });
    });

    // Use the stored pixels to compute the results
    py = 0;
    LOOP(PX_PER_WI_Y, py,
    {
        int y = startY + py;
        int px = 0;
        LOOP(PX_PER_WI_X, px,
        {
            int x = startX + px;
            INTERMEDIATE_TYPE total_sum = 0;
            int sy = 0;
            LOOP(KERNEL_SIZE_Y, sy,
            {
                int sx = 0;
                LOOP(KERNEL_SIZE_X, sx,
                {
                    total_sum = mad(privateKernel[sy][sx], privateData[py + sy][px + sx], total_sum);
                });
            });

            __global TYPE* dstPtr = (__global TYPE*)((__global char*)dstptr + y * dststep + dstoffset + x * sizeof(TYPE)); // Pointer can be out of bounds!
            *dstPtr = CONVERT_TO_TYPE(total_sum);
        });
    });
}
