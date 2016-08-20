/*Daala video codec
Copyright (c) 2015 Daala project contributors.  All rights reserved.

Redistribution and use in source and binary forms, with or without
modification, are permitted provided that the following conditions are met:

- Redistributions of source code must retain the above copyright notice, this
  list of conditions and the following disclaimer.

- Redistributions in binary form must reproduce the above copyright notice,
  this list of conditions and the following disclaimer in the documentation
  and/or other materials provided with the distribution.

THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE
DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR CONTRIBUTORS BE LIABLE
FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL
DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR
SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER
CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY,
OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.*/

/* clang-format off */

#if !defined(_zigzag_H)
# define _zigzag_H (1)

extern const unsigned char OD_ZIGZAG4_DCT_DCT[15][2];
extern const unsigned char OD_ZIGZAG4_ADST_DCT[15][2];
extern const unsigned char OD_ZIGZAG4_DCT_ADST[15][2];
extern const unsigned char OD_ZIGZAG4_ADST_ADST[15][2];

extern const unsigned char OD_ZIGZAG8_DCT_DCT[48][2];
extern const unsigned char OD_ZIGZAG8_ADST_DCT[48][2];
extern const unsigned char OD_ZIGZAG8_DCT_ADST[48][2];
extern const unsigned char OD_ZIGZAG8_ADST_ADST[48][2];

extern const unsigned char OD_ZIGZAG16_DCT_DCT[192][2];
extern const unsigned char OD_ZIGZAG16_ADST_DCT[192][2];
extern const unsigned char OD_ZIGZAG16_DCT_ADST[192][2];
extern const unsigned char OD_ZIGZAG16_ADST_ADST[192][2];

extern const unsigned char OD_ZIGZAG32_DCT_DCT[768][2];

extern const unsigned char OD_ZIGZAG64_DCT_DCT[3072][2];

#endif
