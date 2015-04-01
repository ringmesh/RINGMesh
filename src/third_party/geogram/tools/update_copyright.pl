#!/usr/bin/perl
use strict;
use warnings FATAL => 'all';
use FindBin();
use File::Basename();
use POSIX();
use lib $FindBin::Bin;
use FileUtils();

#############################################################################
# Configuration
#############################################################################

#! Help sections
my $HELP = {
    NAME => <<END,
    $FindBin::Script - updates headers and copyrights in source code
END

    DESCRIPTION => <<END,
    $FindBin::Script searches C and C++ files in the specified file arguments
    and updates the following:
    - leading copyright comment
    - protection macro in header files.

    By default the result of file F is dumped to file F.copyright unless options
    --replace or --no-backup are specified.

    NOTE: Only files with extension (c, h, cpp, hpp, cxx, hxx, C, H) are
    processed.
END
};

#! Option values
my %OPTIONS;

#! Current year
my $YEAR = POSIX::strftime('%Y', localtime());

#! The text of the copyright
my $COPYRIGHT = <<END;
/*
 *      \\V (O |R |P /A |L |I |N |E
 * (C) Bruno Levy, INRIA - ALICE, 2012-$YEAR
 *
 *   Confidential - proprietary software
 */
END

#############################################################################
# Main
#############################################################################

FileUtils::update_files(
    id => 'copyright',
    option_values => \%OPTIONS,
    text_handler => \&process_text,
    help => $HELP,
);

#############################################################################

#!
# @brief Process a single file
# @param[in] file path to the input file
#
sub process_text {
    my($file, $text) = @_;

    # Remove the leading copyright

    $text =~ s{^\s*/\*.+?INRIA.+?\*/\s*}{}s;
    $text =~ s{\s+$}{}s;


    # In header files:
    # - Add the correct copyright
    # - Replace the header protection macro
    # In impl files:
    # - Add the correct copyright
    #
    # NOTE: when the text is read from the standard input, we have no file
    # name information to update the header protection macros, thus we can
    # only update the copyright.

    if( not defined($file) or $file =~ m{\.(?:c|cpp|cxx|C)$} ) {
        $text = <<END;
$COPYRIGHT
$text

END
    }
    else {
        # Remove the existing protection macro (if any)
        # NOTE: the last #endif may have a trailing comment

        if( $text =~ s{^\s*#ifndef (__VORPALINE_\w+?__)\s+#define \1\s+}{}s ) {
            $text =~ s{\s*#endif(?:\s*//\N*)?$}{}s;
        }

        # Add the new protection macro

        my $name = File::Basename::basename($file);
        my $dir = File::Basename::basename(File::Basename::dirname($file));
        my $name_we = $name;
        $name_we =~ s{\.(?:h|hpp|hxx|H)$}{};

        my $protection_macro = uc("__VORPALINE_${dir}_${name_we}__");

        $text = <<END;
$COPYRIGHT
#ifndef $protection_macro
#define $protection_macro

$text

#endif

END
    }

    return $text;
}

