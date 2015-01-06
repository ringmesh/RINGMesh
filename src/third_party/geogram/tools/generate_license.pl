#!/usr/bin/perl
#
# This script generates a license file for Vorpaline
# It can be called either from a build tree or from the
# outside by specifying the path to the desired build tree
#
use strict;
use warnings FATAL => 'all';
use Getopt::Long();
use POSIX();

#!
# \brief Encryption keys
# \details TODO: instead of duplicating the key values,
# we should read the keys from license.h directly 
#
my $VOR_MESSAGE_KEY = "AUATU\0T\$.H\0\\\$\@A\0VUUUUUUUH\0T\$ H\0H+C\0ATUH\0";

#! Output directory where to store generated files
my $OUTPUT_DIR;

# Get command line arguments
Getopt::Long::GetOptions(
    'dir=s' => \$OUTPUT_DIR,
) or exit(1);

#! Duration of the license in days
my $LICENSE_DAYS = shift(@ARGV) || 0;

#! Name of the licensee
my $LICENSE_CLIENT = shift(@ARGV) || "Internal version";

# If specified, make sure that the output dir is a CMake build tree

if( defined($OUTPUT_DIR) ) {
    if( not -f "$OUTPUT_DIR/CMakeCache.txt" ) {
        die "Error: directory $OUTPUT_DIR is not a CMake build tree\n";
    }
}

# Generate the license files
my $BUILD_TIME = time();
generate_license_header();
generate_license_info();
exit(0);


##############################################################################
#!
# \brief Generates license header file license_info.h
#
sub generate_license_header {

    my $license_time = $LICENSE_DAYS*(24*3600);
    my $build_date   = localtime($BUILD_TIME);
    my $max_time     = $BUILD_TIME + $license_time;

    my $text = <<END;
/*
 * License issued to $LICENSE_CLIENT
 * License expires in $LICENSE_DAYS day(s) from $build_date
 */
#ifndef __VORPALINE_LICENSE_INFO__
#define __VORPALINE_LICENSE_INFO__
END

    if( $LICENSE_DAYS == 0 ) {
        $text .= "#define VOR_NO_LICENSE\n";
        $text .= "#define VOR_MAXTIME 0\n";
    } else {
        $text .= "#define VOR_MAXTIME $max_time\n";
        $text .= "#define VOR_CLIENT \"" . encode_string($LICENSE_CLIENT, $VOR_MESSAGE_KEY, 3) . "\"\n";
        $text .= "#define VOR_BUILD_DATE \"" . encode_string($build_date, $VOR_MESSAGE_KEY) . "\"\n";
    }

    $text .= "#endif\n\n";

    if( defined($OUTPUT_DIR) ) {
        save_file("$OUTPUT_DIR/src/lib/vorpalib/license_info.h", $text);
    } else {
        print $text;
    }

}

##############################################################################
#!
# \brief Generates the license information file LICENSE.txt
#
sub generate_license_info {

    my $build_date = localtime($BUILD_TIME);
    my $year = POSIX::strftime('%Y', localtime($BUILD_TIME));

    my $text = <<END;
\\V (O |R |P /A |L |I |N |E
(C) Bruno Levy, INRIA - ALICE, 2012-$year

License issued to $LICENSE_CLIENT
License expires in $LICENSE_DAYS day(s) from $build_date
License is personal and non-transferable
END

    if( defined($OUTPUT_DIR) ) {
	mkdir("$OUTPUT_DIR/doc");
        save_file("$OUTPUT_DIR/doc/LICENSE.txt", $text);
    }
}

##############################################################################
#!
# \brief Encodes the message by XORing message chars with key chars
# \param[in] msg the message to encode
# \param[in] key the encryption key
# \param[in] key_index index of the first key char (default = 0)
#
sub encode_string {
    my($msg, $key, $key_index) = @_;

    if( not defined($key_index) ) {
        $key_index = 0;
    }

    my $msg_length = length($msg);
    my @key_chars = split(//, $key);

    $msg =~ s{.}{ sprintf("\\x%x", ord($&) ^ ord($key_chars[$key_index++])) }ges;
    return sprintf("\\x%x%s", $msg_length, $msg);
}

##############################################################################
#!
# \brief Writes a text to a file
# \param[in] file path to the file to save
# \param[in] text text to write to \p file
#
sub save_file {
    my($file, $text) = @_;

    print "Generating file $file\n";

    if( not open(FILE, '>', $file) ) {
        die "Error: failed to save file $file: $!\n";
    }

    print FILE $text;
    close(FILE);
}


