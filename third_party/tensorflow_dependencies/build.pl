#!/usr/bin/env perl
#
# Copyright (c) 2020, Alliance for Open Media. All rights reserved
#
# This source code is subject to the terms of the BSD 2 Clause License and
# the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
# was not distributed with this source code in the LICENSE file, you can
# obtain it at www.aomedia.org/license/software. If the Alliance for Open
# Media Patent License 1.0 was not distributed with this source code in the
# PATENTS file, you can obtain it at www.aomedia.org/license/patent.
#
###########################################################################
#
# Script to build the static TensorFlow lite library.
#
# Tensorflow's build process generates the static libraries in the same
# directory as the source code. AOM, however, generates binaries/objects/etc.
# in a different directory. This script:
#
# 1.) Copies the TensorFlow code to a temporary directory
# 2.) Copies the necessary dependencies into the temporary directory
# 3.) Compiles it
# 4.) Copies the static library to the AOM build directory
#
# Note that we do not use the download_dependencies.sh script directly, as
# it downloads directly into the source directory.

use strict;
use warnings;
use autodie;
use Archive::Tar;
use Cwd;
use Digest::SHA;
use File::Basename;
use File::Copy;
use File::Copy::Recursive;
use File::Spec::Functions;
use File::Temp;
use LWP::Simple;

# Download the fft2d library, verify the checksum, extract it from the
# tarball, and apply a tweak to it (TF Lite's "download_dependencies.sh"
# performs 3 tweaks, but only 1 of them is relevant).
sub download_fft2d {
    # Download the library.
    my $url = "https://storage.googleapis.com/mirror.tensorflow.org/www.kurims.kyoto-u.ac.jp/~ooura/fft2d.tgz";
    my $downloads_dir = $_[0];
    my $file = catfile($downloads_dir, "fft2d.tgz");
    getstore($url, $file);

    # Check the SHA256 digest.
    my $sha2 = Digest::SHA->new(256);
    $sha2->addfile($file, "b");
    my $digest = $sha2->hexdigest;
    if ($digest ne "ada7e99087c4ed477bfdf11413f2ba8db8a840ba9bbf8ac94f4f3972e2a7cec9") {
        die "Bad checksum for fft2d.tgz: $digest";
    }
    # Extract the code from the tarball.
    my $tar = Archive::Tar->new;
    $tar->read($file);
    my $curr_dir = getcwd;
    chdir $downloads_dir;
    $tar->extract();
    chdir $curr_dir;

    # Apply the tweak.
    my $complex_file = catfile($downloads_dir, "eigen", "Eigen", "src", "Core",
                               "arch", "NEON", "Complex.h");
    open(my $fh, '<', $complex_file);
    my @lines = <$fh>;
    close($fh);
    open($fh, '>', $file);
    foreach my $line (@lines) {
        $line =~ s/static uint64x2_t p2ul_CONJ_XOR = vld1q_u64\( p2ul_conj_XOR_DATA \);/static uint64x2_t p2ul_CONJ_XOR;/;
        print $fh $line;
    }
    close($fh);
}

sub copy_tensorflow_lite_dependencies {
    my $source_dir = $_[0];
    my $output_dir = $_[1];
    my $dependencies_dir = catfile($source_dir, "third_party",
                                   "tensorflow_dependencies");
    opendir(my $dh, $dependencies_dir);
    while (my $file = readdir($dh)) {
        # Ignore files that begin with a period.
        next if ($file =~ m/^\./);
        # Ignore non-directories.
        next unless (-d catfile($dependencies_dir, $file));
        # Copy the directory.
        print "  * $file\n";
        File::Copy::Recursive::dircopy(
          catfile($dependencies_dir, $file),
          catfile($output_dir, $file));
    }
    closedir($dh);
}

# Start of program logic.
if ($#ARGV + 1 != 2) {
    my $prog = basename($0);
    warn("Usage: $prog <source directory> <output directory>\n");
    exit(1);
}
my ($source_dir, $output_dir) = @ARGV;

my $temp_dir = File::Temp::tempdir( CLEANUP => 1 );

print "Copying TensorFlow code to temporary directory for building...\n";
my $tf_dir = catfile($source_dir, "third_party", "tensorflow");
File::Copy::Recursive::dircopy($tf_dir, $temp_dir);

my $downloads_dir = catfile($temp_dir, "tensorflow", "lite", "tools",
                            "make", "downloads");
mkdir $downloads_dir;
print "Copying TensorFlow Lite dependencies to temporary directory...\n";
copy_tensorflow_lite_dependencies($source_dir, $downloads_dir);

print "Downloading fft2d library...\n";
download_fft2d($downloads_dir);

print "Building TensorFlow Lite...\n";
my $build_sh = catfile($temp_dir, "tensorflow", "lite", "tools", "make",
                       "build_lib.sh");
`$build_sh`;
my $liblite = catfile($temp_dir, "tensorflow", "lite", "tools", "make", "gen",
                      "linux_x86_64", "lib", "libtensorflow-lite.a");
print "Copying static library into build directory...\n";
copy($liblite, $output_dir);
