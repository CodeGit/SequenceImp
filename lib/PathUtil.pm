package PathUtil;

=head1 NAME

PathUtil - A support module for checking the PATH variable on *nix and windows OSes

=head1 SYNOPSIS

    use PathUtil qw(findExecutable);
    my $exe = findExecutable("gzip");
    print "Executable: ".$exe."\n";
    my $separator = getDirectorySeparator();
    print "Separator $separator\n";

=head2 Methods

=over 12

=item C<findExecutable>

Returns a string representation of the full path to the string name of the executable supplied as an argument
to the function. The executable will have ".exe" appended in MS Windows.

=item C<getDirectorySeparator>

Returns a string representing the directory separator for the OS. '\' for MS Windows and '/' in Linux and
other unix-based OSes

=back

=head1 LICENSE

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program (Please see the COPYING file for details).
If not, see <http://www.gnu.org/licenses/>.

Please send bug reports to:
kraken@ebi.ac.uk

=head1 AUTHOR

Matloob Qureshi

=cut

use strict;
use warnings;
use Exporter;
use vars qw($VERSION @ISA @EXPORT @EXPORT_OK %EXPORT_TAGS);

$VERSION = 1.00;
@ISA = qw(Exporter);
@EXPORT = ();
@EXPORT_OK = qw(findExecutable getDirectorySeparator);
%EXPORT_TAGS = (DEFAULT => [qw(&findExecutable &getDirectorySeparator)]);

sub findExecutable {
	my $executable = shift @_;
    my $pathStr = $ENV{PATH};
    my @paths;
    my $directorySeparator = getDirectorySeparator();

	if ($^O eq "MSWin32") {
        $executable .= ".exe";
        (@paths) = split(/;/, $pathStr);
    } else {
        (@paths) = split(/:/, $pathStr);
    }

    foreach my $path (@paths) {
    	if ($path !~ /.*?$directorySeparator$/) {
            $path .= $directorySeparator;
        }
        my $exePath = $path . $executable;
        #print STDERR "checking $exePath\n";
        if (-f $exePath and -x $exePath) {
            #$exePath =~ s|\\|\/|g;
            return $exePath;
        }
    }
}

sub getDirectorySeparator {
     if ($^O eq "MSWin32") {
            return "\\";
        } else {
            return "\/";
        }
}

1;