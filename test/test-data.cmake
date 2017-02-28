##
## Copyright (c) 2017, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##

set(AOM_TEST_DATA_LIST "${AOM_ROOT}/test/test-data.sha1")
set(AOM_TEST_DATA_URL "http://downloads.webmproject.org/test_data/libvpx/")
set(AOM_TEST_DATA_PATH "$ENV{LIBAOM_TEST_DATA_PATH}")

if (${AOM_TEST_DATA_PATH} STREQUAL "")
  message(WARNING "Writing test data to ${AOM_CONFIG_DIR}, set "
          "$LIBAOM_TEST_DATA_PATH in your environment to avoid this warning.")
  set(AOM_TEST_DATA_PATH "${AOM_CONFIG_DIR}")
endif ()

if (NOT EXISTS "${AOM_TEST_DATA_PATH}")
  file(MAKE_DIRECTORY "${AOM_TEST_DATA_PATH}")
endif ()

# Read test-data.sha1 into $files_and_checksums. $files_and_checksums becomes a
# list with an entry for each line from $AOM_TEST_DATA_LIST.
file(STRINGS "${AOM_TEST_DATA_LIST}" files_and_checksums)

# Iterate over the list of lines and split it into $checksums and $filenames.
foreach (line ${files_and_checksums})
  string(FIND "${line}" " *" delim_pos)

  math(EXPR filename_pos "${delim_pos} + 2")
  string(SUBSTRING "${line}" 0 ${delim_pos} checksum)
  string(SUBSTRING "${line}" ${filename_pos} -1 filename)

  set(checksums ${checksums} ${checksum})
  set(filenames ${filenames} ${filename})
endforeach ()

message("----> checksums=${checksums}")
message("----> filenames=${filenames}")


list(LENGTH filenames num_files)
math(EXPR num_files "${num_files} - 1")

foreach (file_num RANGE ${num_files})
  list(GET filenames ${file_num} filename)
  list(GET checksums ${file_num} checksum)
  message("----> ${file_num} ${filename} ${checksum}")
  set(filename "${AOM_TEST_DATA_PATH}/${filename}")
  if (EXISTS "${filename}")
    file(SHA1 "${filename}" file_sha1)
  endif ()
  if (NOT "${file_sha1}" STREQUAL "${checksum}")
    # checksum mismatch or file doesn't exist; download it.
    message("----> checksum mismatch for ${filename}.")
    message("----> expected |${checksum}| got |${sha1sum}|")
  endif ()
endforeach ()
