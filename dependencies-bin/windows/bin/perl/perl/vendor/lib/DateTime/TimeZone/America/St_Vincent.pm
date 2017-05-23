# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/6Pwc8w6J1M/northamerica.  Olson data version 2013d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::St_Vincent;
{
  $DateTime::TimeZone::America::St_Vincent::VERSION = '1.60';
}
BEGIN {
  $DateTime::TimeZone::America::St_Vincent::AUTHORITY = 'cpan:DROLSKY';
}

use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::St_Vincent::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59611176296, #      utc_end 1890-01-01 04:04:56 (Wed)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59611161600, #    local_end 1890-01-01 00:00:00 (Wed)
-14696,
0,
'LMT',
    ],
    [
59611176296, #    utc_start 1890-01-01 04:04:56 (Wed)
60305313896, #      utc_end 1912-01-01 04:04:56 (Mon)
59611161600, #  local_start 1890-01-01 00:00:00 (Wed)
60305299200, #    local_end 1912-01-01 00:00:00 (Mon)
-14696,
0,
'KMT',
    ],
    [
60305313896, #    utc_start 1912-01-01 04:04:56 (Mon)
DateTime::TimeZone::INFINITY, #      utc_end
60305299496, #  local_start 1912-01-01 00:04:56 (Mon)
DateTime::TimeZone::INFINITY, #    local_end
-14400,
0,
'AST',
    ],
];

sub olson_version { '2013d' }

sub has_dst_changes { 0 }

sub _max_year { 2023 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}



1;

