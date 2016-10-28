##
## Copyright (c) 2016, Alliance for Open Media. All rights reserved
##
## This source code is subject to the terms of the BSD 2 Clause License and
## the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
## was not distributed with this source code in the LICENSE file, you can
## obtain it at www.aomedia.org/license/software. If the Alliance for Open
## Media Patent License 1.0 was not distributed with this source code in the
## PATENTS file, you can obtain it at www.aomedia.org/license/patent.
##
#!/usr/bin/env perl

use strict;
use warnings;
use 5.010;
use Getopt::Long;

my $changelog_filename;
my $version_filename;
GetOptions('changelog_filename=s' => \$changelog_filename,
           'version_filename=s' => \$version_filename);

print("changelog_filename=$changelog_filename\n");
print("version_filename=$version_filename\n");

open(my $changelog_file, '<', $changelog_filename) or
  die("Unable to open CHANGELOG @ $changelog_filename: $!.");

my $version_string;
while (my $line = <$changelog_file>) {
  my @split_line = split(" ", $line, 3);
  $version_string = "" . $split_line[1];
  last if substr($version_string, 0, 1) eq "v";
}
close($changelog_file);

my @version_components = split('\.', $version_string, 4);
my $version_major = substr($version_components[0], 1);
my $version_minor = $version_components[1];
my $version_patch = $version_components[2];

my $version_extra = "";
if (@version_components > 3) {
  $version_extra = $version_components[3];
}

open(my $version_file, '>', $version_filename) or
  die("Cannot open $version_filename: $!");

my $version_packed = "((VERSION_MAJOR<<16)|(VERSION_MINOR<<8)|(VERSION_PATCH))";
my $lic_block = << 'EOF';
/*
 * Copyright (c) 2016, Alliance for Open Media. All rights reserved
 *
 * This source code is subject to the terms of the BSD 2 Clause License and
 * the Alliance for Open Media Patent License 1.0. If the BSD 2 Clause License
 * was not distributed with this source code in the LICENSE file, you can
 * obtain it at www.aomedia.org/license/software. If the Alliance for Open
 * Media Patent License 1.0 was not distributed with this source code in the
 * PATENTS file, you can obtain it at www.aomedia.org/license/patent.
 */
EOF

$version_file->print($lic_block);
$version_file->print("\#define VERSION_MAJOR $version_major\n");
$version_file->print("\#define VERSION_MINOR $version_minor\n");
$version_file->print("\#define VERSION_PATCH $version_patch\n");
$version_file->print("\#define VERSION_EXTRA \"$version_extra\"\n");
$version_file->print("\#define VERSION_PACKED $version_packed\n");
$version_file->print("\#define VERSION_STRING_NOSP \"$version_string\"\n");
$version_file->print("\#define VERSION_STRING \" $version_string\"\n");
close($version_file);
