#AV1 Codec Library

##Building the library and applications

###Prerequisites
 1. [CMake](https://cmake.org) version 3.5 or higher.
 2. [Git](https://git-scm.com/).
 3. [Perl](https://www.perl.org/).
 4. For x86 targets, [yasm](http://yasm.tortall.net/), which is preferred, or a recent version of [nasm](http://www.nasm.us/).
 5. Building the documentation requires [doxygen](http://doxygen.org).

###Basic build
CMake replaces the configure step typical of many projects. Running CMake will produce configuration and build files for the currently selected CMake generator. For most system the default generator is Unix Makefiles. The basic form of a makefile build is the following:

    $ cmake path/to/aom
    $ make

The above will generate a makefile build that produces the AV1 library and applications for the current host system after the make step completes successfully. The compiler chosen varies by host platform, but a general rule applies: On systems where cc and c++ are present in $PATH at the time CMake is run the generated build will use cc and c++ by default.

###Configuration options
The AV1 codec library has a great many configuration options. These come in two varieties:

 1. Build system configuration options. These have the form ENABLE_FEATURE.
 2. AV1 codec configuration options. These have the form CONFIG_FEATURE.

Both types of option are set at the time CMake is run. The following example enables ccache and disables high bit depth:

    $ cmake path/to/aom -DENABLE_CCACHE=1 -DCONFIG_HIGHBITDEPTH=0
    $ make

The available configuration options are too numerous to list here. Build system configuration options can be found at the top of the CMakeLists.txt file found in the root of the AV1 repository, and AV1 codec configuration options can currently be found in the file build/cmake/aom\_config_defaults.cmake.

###Cross compiling
For the purposes of building the AV1 codec and applications and relative to the scope of this guide, all builds for architectures differing from the native host architecture will be considered cross compiles. The AV1 CMake build handles cross compiling via the use of toolchain files included in the AV1 repository. The toolchain files available at the time of this writing are:

 - arm64-ios.cmake
 - arm64-linux-gcc.cmake
 - armv7-ios.cmake
 - armv7-linux-gcc.cmake
 - armv7s-ios.cmake
 - ios-simulator-common.cmake
 - mips32-linux-gcc.cmake
 - mips64-linux-gcc.cmake
 - x86-ios-simulator.cmake
 - x86-linux.cmake
 - x86-macos.cmake
 - x86_64-ios-simulator.cmake

The following example demonstrates use of the x86-macos.cmake toolchain file on a x86_64 MacOS host:

    $ cmake path/to/aom \
    -DCMAKE_TOOLCHAIN_FILE=path/to/aom/build/cmake/toolchains/x86-macos.cmake
    $ make

To build for an unlisted target creation of a new toolchain file is the best solution. The existing toolchain files can be used a starting point for a new toolchain file since each one exposes the basic requirements for toolchain files as used in the AV1 codec build. 

As a temporary work around an unoptimized AV1 configuration that builds only C and C++ sources can be produced using the following commands:

    $ cmake path/to/aom -DAOM_TARGET_CPU=generic
    $ make

In addition to the above it's important to note that the toolchain files suffixed with gcc behave differently than the others. These toolchain files attempt to obey the $CROSS environment variable.

##Testing the AV1 codec

Currently there are two types of tests in the AV1 codec repository:

 1. Unit tests.
 2. Example tests.

The unit tests can be run at build time:

    # Before running the make command the LIBAOM_TEST_DATA_PATH environment
    # variable should be set to avoid downloading the test files to the
    # cmake build configuration directory.
    $ cmake path/to/aom
    # Note: The AV1 CMake build creates many test targets. Running make
    # with multiple jobs will speed up the test run significantly.
    $ make runtests

The example tests require a bash shell and can be run in the following manner:

    # See the note above about LIBAOM_TEST_DATA_PATH above.
    $ cmake path/to/aom
    $ make
    # It's best to build the testdata target using many make jobs.
    # Running it like this will verify and download (if necessary)
    # one at a time, which takes a while.
    $ make testdata
