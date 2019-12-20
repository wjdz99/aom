# AV2 Tensorflow-Lite Integration

Tensorflow-lite is a submodule of the experimental branch of libaom. In order
to access it, run:

git submodule update --init

The TF-lite library has two characteristics that require workaround:

  1. It compiles the TF-lite library and intermediate files in the source
     directory. To prevent the libaom source directory from pollution, the
     script `build.sh` will make a copy of the library, build it, and copy over
     the compiled output into a destination directory.
  2. TF-lite requires the execution of
     `tensorflow/lite/tools/make/download_dependencies.sh` to download
     dependencies. This is an expensive operation, and it does not appear to
     validate the downloaded code. To ensure consistent releases, the downloaded
     code is bundled into downloads.tar.xz, which is copied over by the
     `build.sh` script.
