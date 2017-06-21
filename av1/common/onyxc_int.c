#include <assert.h>
#include "av1/common/onyxc_int.h"

#if CONFIG_HASH_ME
void crc_calculator_light_init(crc_calculator_light *p_crc_calculator_light,
                               uint32_t bits, uint32_t truncPoly) {
  p_crc_calculator_light->m_remainder = 0;
  p_crc_calculator_light->m_bits = bits;
  p_crc_calculator_light->m_truncPoly = truncPoly;
  p_crc_calculator_light->m_finalResultMask = (1 << bits) - 1;
  crc_calculator_light_init_table(p_crc_calculator_light);
}

void crc_calculator_light_process_data(
    crc_calculator_light *p_crc_calculator_light, uint8_t *pData,
    uint32_t dataLength) {
  for (uint32_t i = 0; i < dataLength; i++) {
    uint8_t index = (p_crc_calculator_light->m_remainder >>
                     (p_crc_calculator_light->m_bits - 8)) ^
                    pData[i];
    p_crc_calculator_light->m_remainder <<= 8;
    p_crc_calculator_light->m_remainder ^=
        p_crc_calculator_light->m_table[index];
  }
}

void crc_calculator_light_reset(crc_calculator_light *p_crc_calculator_light) {
  p_crc_calculator_light->m_remainder = 0;
}

uint32_t crc_calculator_light_get_crc(
    crc_calculator_light *p_crc_calculator_light) {
  return p_crc_calculator_light->m_remainder &
         p_crc_calculator_light->m_finalResultMask;
}

void crc_calculator_light_init_table(
    crc_calculator_light *p_crc_calculator_light) {
  const uint32_t highBit = 1 << (p_crc_calculator_light->m_bits - 1);
  const uint32_t ByteHighBit = 1 << (8 - 1);

  for (uint32_t value = 0; value < 256; value++) {
    uint32_t remainder = 0;
    for (uint8_t mask = ByteHighBit; mask != 0; mask >>= 1) {
      if (value & mask) {
        remainder ^= highBit;
      }

      if (remainder & highBit) {
        remainder <<= 1;
        remainder ^= p_crc_calculator_light->m_truncPoly;
      } else {
        remainder <<= 1;
      }
    }
    p_crc_calculator_light->m_table[value] = remainder;
  }
}

static const int m_CRCBits = 16;
static const int m_blockSizeBits = 2;
static int m_blockSizeToIndex[65][65];

static crc_calculator_light m_crcCalculator1;
static crc_calculator_light m_crcCalculator2;
static int g_crcInitialized = 0;

void hash_table_init(hash_table *p_hash_table) {
  if (g_crcInitialized == 0) {
    crc_calculator_light_init(&m_crcCalculator1, 24, 0x5D6DCB);
    crc_calculator_light_init(&m_crcCalculator2, 24, 0x864CFB);
    hash_init_block_size_to_index();
    g_crcInitialized = 1;
  }
  p_hash_table->p_lookup_table = NULL;
}

void hash_table_destroy(hash_table *p_hash_table) {
  hash_table_clearAll(p_hash_table);
  if (p_hash_table->p_lookup_table != NULL) {
    free(p_hash_table->p_lookup_table);
    p_hash_table->p_lookup_table = NULL;
  }
}

void hash_table_create(hash_table *p_hash_table) {
  if (p_hash_table->p_lookup_table != NULL) {
    hash_table_clearAll(p_hash_table);
    return;
  }
  int maxAddr = 1 << (m_CRCBits + m_blockSizeBits);
  p_hash_table->p_lookup_table = (Vector **)malloc(sizeof(Vector *) * maxAddr);
  memset(p_hash_table->p_lookup_table, 0, sizeof(Vector *) * maxAddr);
}

void hash_table_clearAll(hash_table *p_hash_table) {
  if (p_hash_table->p_lookup_table == NULL) {
    return;
  }
  int maxAddr = 1 << (m_CRCBits + m_blockSizeBits);
  for (int i = 0; i < maxAddr; i++) {
    if (p_hash_table->p_lookup_table[i] != NULL) {
      vector_destroy(p_hash_table->p_lookup_table[i]);
      free(p_hash_table->p_lookup_table[i]);
      p_hash_table->p_lookup_table[i] = NULL;
    }
  }
}

void hash_table_add_to_table(hash_table *p_hash_table, uint32_t hashValue,
                             block_hash *blockHash) {
  if (p_hash_table->p_lookup_table[hashValue] == NULL) {
    p_hash_table->p_lookup_table[hashValue] = malloc(sizeof(Vector));
    vector_setup(p_hash_table->p_lookup_table[hashValue], 10,
                 sizeof(block_hash));
    vector_push_back(p_hash_table->p_lookup_table[hashValue], blockHash);
  } else {
    vector_push_back(p_hash_table->p_lookup_table[hashValue], blockHash);
  }
}

int32_t hash_table_count(hash_table *p_hash_table, uint32_t hashValue) {
  if (p_hash_table->p_lookup_table[hashValue] == NULL) {
    return 0;
  } else {
    return (int32_t)(p_hash_table->p_lookup_table[hashValue]->size);
  }
}

Iterator hash_get_first_iterator(hash_table *p_hash_table, uint32_t hashValue) {
  assert(hash_table_count(p_hash_table, hashValue) > 0);
  return vector_begin(p_hash_table->p_lookup_table[hashValue]);
}

int32_t has_exact_match(hash_table *p_hash_table, uint32_t hashValue1,
                        uint32_t hash_value2) {
  if (p_hash_table->p_lookup_table[hashValue1] == NULL) {
    return 0;
  }
  Iterator iterator = vector_begin(p_hash_table->p_lookup_table[hashValue1]);
  Iterator last = vector_end(p_hash_table->p_lookup_table[hashValue1]);
  for (; !iterator_equals(&iterator, &last); iterator_increment(&iterator)) {
    if ((*(block_hash *)iterator_get(&iterator)).hash_value2 == hash_value2) {
      return 1;
    }
  }
  return 0;
}

void generate_block_2x2_hash_value(YV12_BUFFER_CONFIG *picture,
                                   uint32_t *picBlockHash[2],
                                   int8_t *picBlockSameInfo[3]) {
  const int width = 2;
  const int height = 2;
  int xEnd = picture->y_crop_width - width + 1;
  int yEnd = picture->y_crop_height - height + 1;

  int length = width * 2;
  uint8_t *p = malloc(sizeof(uint8_t) * length);

  int pos = 0;
  for (int yPos = 0; yPos < yEnd; yPos++) {
    for (int xPos = 0; xPos < xEnd; xPos++) {
      get_pixels_in_1D_char_array_by_block_2x2(
          picture->y_buffer + yPos * picture->y_stride + xPos,
          picture->y_stride, p);
      picBlockSameInfo[0][pos] = is_block_2x2_row_same_value(p);
      picBlockSameInfo[1][pos] = is_block_2x2_col_same_value(p);

      picBlockHash[0][pos] = get_crc_value1(p, length * sizeof(uint8_t));
      picBlockHash[1][pos] = get_crc_value2(p, length * sizeof(uint8_t));

      pos++;
    }
    pos += width - 1;
  }

  free(p);
}

void generate_block_hash_value(YV12_BUFFER_CONFIG *picture, int width,
                               int height, uint32_t *srcPicBlockHash[2],
                               uint32_t *dstPicBlockHash[2],
                               int8_t *srcPicBlockSameInfo[3],
                               int8_t *dstPicBlockSameInfo[3]) {
  int picWidth = picture->y_crop_width;
  int xEnd = picture->y_crop_width - width + 1;
  int yEnd = picture->y_crop_height - height + 1;

  int srcWidth = width >> 1;
  int quadWidth = width >> 2;
  int srcHeight = height >> 1;
  int quadHeight = height >> 2;

  int length = 4 * sizeof(uint32_t);

  uint32_t *p = malloc(sizeof(uint32_t) * 4);
  int pos = 0;
  for (int yPos = 0; yPos < yEnd; yPos++) {
    for (int xPos = 0; xPos < xEnd; xPos++) {
      p[0] = srcPicBlockHash[0][pos];
      p[1] = srcPicBlockHash[0][pos + srcWidth];
      p[2] = srcPicBlockHash[0][pos + srcHeight * picWidth];
      p[3] = srcPicBlockHash[0][pos + srcHeight * picWidth + srcWidth];
      dstPicBlockHash[0][pos] = get_crc_value1((uint8_t *)p, length);

      p[0] = srcPicBlockHash[1][pos];
      p[1] = srcPicBlockHash[1][pos + srcWidth];
      p[2] = srcPicBlockHash[1][pos + srcHeight * picWidth];
      p[3] = srcPicBlockHash[1][pos + srcHeight * picWidth + srcWidth];
      dstPicBlockHash[1][pos] = get_crc_value2((uint8_t *)p, length);

      dstPicBlockSameInfo[0][pos] =
          srcPicBlockSameInfo[0][pos] &&
          srcPicBlockSameInfo[0][pos + quadWidth] &&
          srcPicBlockSameInfo[0][pos + srcWidth] &&
          srcPicBlockSameInfo[0][pos + srcHeight * picWidth] &&
          srcPicBlockSameInfo[0][pos + srcHeight * picWidth + quadWidth] &&
          srcPicBlockSameInfo[0][pos + srcHeight * picWidth + srcWidth];

      dstPicBlockSameInfo[1][pos] =
          srcPicBlockSameInfo[1][pos] &&
          srcPicBlockSameInfo[1][pos + srcWidth] &&
          srcPicBlockSameInfo[1][pos + quadHeight * picWidth] &&
          srcPicBlockSameInfo[1][pos + quadHeight * picWidth + srcWidth] &&
          srcPicBlockSameInfo[1][pos + srcHeight * picWidth] &&
          srcPicBlockSameInfo[1][pos + srcHeight * picWidth + srcWidth];
      pos++;
    }
    pos += width - 1;
  }

  if (width >= 8) {
    int widthMinus1 = width - 1;
    int heightMinus1 = height - 1;
    pos = 0;
    for (int yPos = 0; yPos < yEnd; yPos++) {
      for (int xPos = 0; xPos < xEnd; xPos++) {
        dstPicBlockSameInfo[2][pos] =
            (!dstPicBlockSameInfo[0][pos] && !dstPicBlockSameInfo[1][pos]) ||
            (((xPos & widthMinus1) == 0) && ((yPos & heightMinus1) == 0));
        pos++;
      }
      pos += width - 1;
    }
  }

  free(p);
}

void add_to_hash_map_by_row_with_precal_data(hash_table *p_hash_table,
                                             uint32_t *picHash[2],
                                             int8_t *picIsSame, int picWidth,
                                             int picHeight, int width,
                                             int height) {
  int xEnd = picWidth - width + 1;
  int yEnd = picHeight - height + 1;

  int8_t *srcIsAdded = picIsSame;
  uint32_t *srcHash[2] = { picHash[0], picHash[1] };

  int addValue = m_blockSizeToIndex[width][height];
  assert(addValue >= 0);
  addValue <<= m_CRCBits;
  int crcMask = 1 << m_CRCBits;
  crcMask -= 1;

  for (int xPos = 0; xPos < xEnd; xPos++) {
    for (int yPos = 0; yPos < yEnd; yPos++) {
      int pos = yPos * picWidth + xPos;
      // valid data
      if (srcIsAdded[pos]) {
        block_hash blockHash;
        blockHash.x = xPos;
        blockHash.y = yPos;

        uint32_t hashValue1 = (srcHash[0][pos] & crcMask) + addValue;
        blockHash.hash_value2 = srcHash[1][pos];

        hash_table_add_to_table(p_hash_table, hashValue1, &blockHash);
      }
    }
  }
}

uint32_t get_crc_value1(uint8_t *p, int length) {
  crc_calculator_light_reset(&m_crcCalculator1);
  crc_calculator_light_process_data(&m_crcCalculator1, p, length);
  return crc_calculator_light_get_crc(&m_crcCalculator1);
}

uint32_t get_crc_value2(uint8_t *p, int length) {
  crc_calculator_light_reset(&m_crcCalculator2);
  crc_calculator_light_process_data(&m_crcCalculator2, p, length);
  return crc_calculator_light_get_crc(&m_crcCalculator2);
}

void get_pixels_in_1D_char_array_by_block_2x2(uint8_t *ySrc, int stride,
                                              uint8_t *pPixelsIn1D) {
  uint8_t *pPel = ySrc;
  int index = 0;
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      pPixelsIn1D[index++] = pPel[j];
    }
    pPel += stride;
  }
}

int is_block_2x2_row_same_value(uint8_t *p) {
  if (p[0] != p[1] || p[2] != p[3]) {
    return 0;
  }

  return 1;
}

int is_block_2x2_col_same_value(uint8_t *p) {
  if ((p[0] != p[2]) || (p[1] != p[3])) {
    return 0;
  }

  return 1;
}

int hash_is_horizontal_perfect(YV12_BUFFER_CONFIG *picture, int width,
                               int height, int xStart, int yStart) {
  int stride = picture->y_stride;
  uint8_t *p = picture->y_buffer;
  p += (yStart)*stride + (xStart);

  for (int i = 0; i < height; i++) {
    for (int j = 1; j < width; j++) {
      if (p[j] != p[0]) {
        return 0;
      }
    }
    p += stride;
  }

  return 1;
}

int hash_is_vertical_perfect(YV12_BUFFER_CONFIG *picture, int width, int height,
                             int xStart, int yStart) {
  int stride = picture->y_stride;
  uint8_t *p = picture->y_buffer;
  p += (yStart)*stride + (xStart);

  for (int i = 0; i < width; i++) {
    for (int j = 1; j < height; j++) {
      if (p[j * stride + i] != p[i]) {
        return 0;
      }
    }
  }

  return 1;
}

int get_block_hash_value(uint8_t *ySrc, int stride, int width, int height,
                         uint32_t *hashValue1, uint32_t *hash_value2) {
  int addValue = m_blockSizeToIndex[width][height];
  assert(addValue >= 0);
  addValue <<= m_CRCBits;
  int crcMask = 1 << m_CRCBits;
  crcMask -= 1;
  int length = 4;

  uint8_t *p = malloc(sizeof(uint8_t) * length);
  uint32_t *toHash = malloc(sizeof(uint32_t) * 4);

  int block2x2Num = (width * height) >> 2;

  uint32_t *hashValueBuffer[2][2];
  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      hashValueBuffer[i][j] = malloc(sizeof(uint32_t) * block2x2Num);
    }
  }

  // 2x2 subblock hash values in current CU
  int subBlockInWidth = (width >> 1);
  for (int yPos = 0; yPos < height; yPos += 2) {
    for (int xPos = 0; xPos < width; xPos += 2) {
      int pos = (yPos >> 1) * subBlockInWidth + (xPos >> 1);
      get_pixels_in_1D_char_array_by_block_2x2(ySrc + yPos * stride + xPos,
                                               stride, p);

      hashValueBuffer[0][0][pos] = get_crc_value1(p, length * sizeof(uint8_t));
      hashValueBuffer[1][0][pos] = get_crc_value2(p, length * sizeof(uint8_t));
    }
  }

  int srcSubBlockInWidth = subBlockInWidth;
  subBlockInWidth >>= 1;
  length = 4 * sizeof(uint32_t);

  int srcIdx = 1;
  int dstIdx = 0;

  // 4x4 subblock hash values to current block hash values
  for (int subWidth = 4; subWidth <= width; subWidth *= 2) {
    srcIdx = 1 - srcIdx;
    dstIdx = 1 - dstIdx;

    int dstPos = 0;
    for (int yPos = 0; yPos < subBlockInWidth; yPos++) {
      for (int xPos = 0; xPos < subBlockInWidth; xPos++) {
        int srcPos = (yPos << 1) * srcSubBlockInWidth + (xPos << 1);

        toHash[0] = hashValueBuffer[0][srcIdx][srcPos];
        toHash[1] = hashValueBuffer[0][srcIdx][srcPos + 1];
        toHash[2] = hashValueBuffer[0][srcIdx][srcPos + srcSubBlockInWidth];
        toHash[3] = hashValueBuffer[0][srcIdx][srcPos + srcSubBlockInWidth + 1];

        hashValueBuffer[0][dstIdx][dstPos] =
            get_crc_value1((uint8_t *)toHash, length);

        toHash[0] = hashValueBuffer[1][srcIdx][srcPos];
        toHash[1] = hashValueBuffer[1][srcIdx][srcPos + 1];
        toHash[2] = hashValueBuffer[1][srcIdx][srcPos + srcSubBlockInWidth];
        toHash[3] = hashValueBuffer[1][srcIdx][srcPos + srcSubBlockInWidth + 1];
        hashValueBuffer[1][dstIdx][dstPos] =
            get_crc_value2((uint8_t *)toHash, length);
        dstPos++;
      }
    }

    srcSubBlockInWidth = subBlockInWidth;
    subBlockInWidth >>= 1;
  }

  *hashValue1 = (hashValueBuffer[0][dstIdx][0] & crcMask) + addValue;
  *hash_value2 = hashValueBuffer[1][dstIdx][0];

  free(toHash);

  for (int i = 0; i < 2; i++) {
    for (int j = 0; j < 2; j++) {
      free(hashValueBuffer[i][j]);
    }
  }

  free(p);

  return 1;
}

void hash_init_block_size_to_index() {
  for (int i = 0; i < 65; i++) {
    for (int j = 0; j < 65; j++) {
      m_blockSizeToIndex[i][j] = -1;
    }
  }

  m_blockSizeToIndex[8][8] = 0;
  m_blockSizeToIndex[16][16] = 1;
  m_blockSizeToIndex[32][32] = 2;
  m_blockSizeToIndex[64][64] = 3;
}

#endif
