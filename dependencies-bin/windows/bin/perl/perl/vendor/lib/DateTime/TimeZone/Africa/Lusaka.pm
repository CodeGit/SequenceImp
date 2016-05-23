# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/6Pwc8w6J1M/africa.  Olson data version 2013d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Africa::Lusaka;
{
  $DateTime::TimeZone::Africa::Lusaka::VERSION = '1.60';
}
BEGIN {
  $DateTime::TimeZone::Africa::Lusaka::AUTHORITY = 'cpan:DROLSKY';
}

use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Africa::Lusaka::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
60026393212, #      utc_end 1903-02-28 22:06:52 (Sat)
DateTime::TimeZone::NEG_INFINITY, #  local_start
60026400000, #    local_end 1903-03-01 00:00:00 (Sun)
6788,
0,
'LMT',
    ],
    [
60026393212, #    utc_start 1903-02-28 22:06:52 (Sat)
DateTime::TimeZone::INFINITY, #      utc_end
60026400412, #  local_start 1903-03-01 00:06:52 (Sun)
DateTime::TimeZone::INFINITY, #    local_end
7200,
0,
'CAT',
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

