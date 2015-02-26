# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/6Pwc8w6J1M/southamerica.  Olson data version 2013d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::Santiago;
{
  $DateTime::TimeZone::America::Santiago::VERSION = '1.60';
}
BEGIN {
  $DateTime::TimeZone::America::Santiago::AUTHORITY = 'cpan:DROLSKY';
}

use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::Santiago::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59611178566, #      utc_end 1890-01-01 04:42:46 (Wed)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59611161600, #    local_end 1890-01-01 00:00:00 (Wed)
-16966,
0,
'LMT',
    ],
    [
59611178566, #    utc_start 1890-01-01 04:42:46 (Wed)
60242244166, #      utc_end 1910-01-01 04:42:46 (Sat)
59611161600, #  local_start 1890-01-01 00:00:00 (Wed)
60242227200, #    local_end 1910-01-01 00:00:00 (Sat)
-16966,
0,
'SMT',
    ],
    [
60242244166, #    utc_start 1910-01-01 04:42:46 (Sat)
60447272400, #      utc_end 1916-07-01 05:00:00 (Sat)
60242226166, #  local_start 1909-12-31 23:42:46 (Fri)
60447254400, #    local_end 1916-07-01 00:00:00 (Sat)
-18000,
0,
'CLT',
    ],
    [
60447272400, #    utc_start 1916-07-01 05:00:00 (Sat)
60515700166, #      utc_end 1918-09-01 04:42:46 (Sun)
60447255434, #  local_start 1916-07-01 00:17:14 (Sat)
60515683200, #    local_end 1918-09-01 00:00:00 (Sun)
-16966,
0,
'SMT',
    ],
    [
60515700166, #    utc_start 1918-09-01 04:42:46 (Sun)
60541876800, #      utc_end 1919-07-01 04:00:00 (Tue)
60515685766, #  local_start 1918-09-01 00:42:46 (Sun)
60541862400, #    local_end 1919-07-01 00:00:00 (Tue)
-14400,
0,
'CLT',
    ],
    [
60541876800, #    utc_start 1919-07-01 04:00:00 (Tue)
60799696966, #      utc_end 1927-09-01 04:42:46 (Thu)
60541859834, #  local_start 1919-06-30 23:17:14 (Mon)
60799680000, #    local_end 1927-09-01 00:00:00 (Thu)
-16966,
0,
'SMT',
    ],
    [
60799696966, #    utc_start 1927-09-01 04:42:46 (Thu)
60818097600, #      utc_end 1928-04-01 04:00:00 (Sun)
60799682566, #  local_start 1927-09-01 00:42:46 (Thu)
60818083200, #    local_end 1928-04-01 00:00:00 (Sun)
-14400,
1,
'CLST',
    ],
    [
60818097600, #    utc_start 1928-04-01 04:00:00 (Sun)
60831320400, #      utc_end 1928-09-01 05:00:00 (Sat)
60818079600, #  local_start 1928-03-31 23:00:00 (Sat)
60831302400, #    local_end 1928-09-01 00:00:00 (Sat)
-18000,
0,
'CLT',
    ],
    [
60831320400, #    utc_start 1928-09-01 05:00:00 (Sat)
60849633600, #      utc_end 1929-04-01 04:00:00 (Mon)
60831306000, #  local_start 1928-09-01 01:00:00 (Sat)
60849619200, #    local_end 1929-04-01 00:00:00 (Mon)
-14400,
1,
'CLST',
    ],
    [
60849633600, #    utc_start 1929-04-01 04:00:00 (Mon)
60862856400, #      utc_end 1929-09-01 05:00:00 (Sun)
60849615600, #  local_start 1929-03-31 23:00:00 (Sun)
60862838400, #    local_end 1929-09-01 00:00:00 (Sun)
-18000,
0,
'CLT',
    ],
    [
60862856400, #    utc_start 1929-09-01 05:00:00 (Sun)
60881169600, #      utc_end 1930-04-01 04:00:00 (Tue)
60862842000, #  local_start 1929-09-01 01:00:00 (Sun)
60881155200, #    local_end 1930-04-01 00:00:00 (Tue)
-14400,
1,
'CLST',
    ],
    [
60881169600, #    utc_start 1930-04-01 04:00:00 (Tue)
60894392400, #      utc_end 1930-09-01 05:00:00 (Mon)
60881151600, #  local_start 1930-03-31 23:00:00 (Mon)
60894374400, #    local_end 1930-09-01 00:00:00 (Mon)
-18000,
0,
'CLT',
    ],
    [
60894392400, #    utc_start 1930-09-01 05:00:00 (Mon)
60912705600, #      utc_end 1931-04-01 04:00:00 (Wed)
60894378000, #  local_start 1930-09-01 01:00:00 (Mon)
60912691200, #    local_end 1931-04-01 00:00:00 (Wed)
-14400,
1,
'CLST',
    ],
    [
60912705600, #    utc_start 1931-04-01 04:00:00 (Wed)
60925928400, #      utc_end 1931-09-01 05:00:00 (Tue)
60912687600, #  local_start 1931-03-31 23:00:00 (Tue)
60925910400, #    local_end 1931-09-01 00:00:00 (Tue)
-18000,
0,
'CLT',
    ],
    [
60925928400, #    utc_start 1931-09-01 05:00:00 (Tue)
60944328000, #      utc_end 1932-04-01 04:00:00 (Fri)
60925914000, #  local_start 1931-09-01 01:00:00 (Tue)
60944313600, #    local_end 1932-04-01 00:00:00 (Fri)
-14400,
1,
'CLST',
    ],
    [
60944328000, #    utc_start 1932-04-01 04:00:00 (Fri)
60957550800, #      utc_end 1932-09-01 05:00:00 (Thu)
60944310000, #  local_start 1932-03-31 23:00:00 (Thu)
60957532800, #    local_end 1932-09-01 00:00:00 (Thu)
-18000,
0,
'CLT',
    ],
    [
60957550800, #    utc_start 1932-09-01 05:00:00 (Thu)
61265131200, #      utc_end 1942-06-01 04:00:00 (Mon)
60957536400, #  local_start 1932-09-01 01:00:00 (Thu)
61265116800, #    local_end 1942-06-01 00:00:00 (Mon)
-14400,
1,
'CLST',
    ],
    [
61265131200, #    utc_start 1942-06-01 04:00:00 (Mon)
61270405200, #      utc_end 1942-08-01 05:00:00 (Sat)
61265113200, #  local_start 1942-05-31 23:00:00 (Sun)
61270387200, #    local_end 1942-08-01 00:00:00 (Sat)
-18000,
0,
'CLT',
    ],
    [
61270405200, #    utc_start 1942-08-01 05:00:00 (Sat)
61395163200, #      utc_end 1946-07-15 04:00:00 (Mon)
61270390800, #  local_start 1942-08-01 01:00:00 (Sat)
61395148800, #    local_end 1946-07-15 00:00:00 (Mon)
-14400,
1,
'CLST',
    ],
    [
61395163200, #    utc_start 1946-07-15 04:00:00 (Mon)
61399306800, #      utc_end 1946-09-01 03:00:00 (Sun)
61395148800, #  local_start 1946-07-15 00:00:00 (Mon)
61399292400, #    local_end 1946-08-31 23:00:00 (Sat)
-14400,
1,
'CLST',
    ],
    [
61399306800, #    utc_start 1946-09-01 03:00:00 (Sun)
61417627200, #      utc_end 1947-04-01 04:00:00 (Tue)
61399288800, #  local_start 1946-08-31 22:00:00 (Sat)
61417609200, #    local_end 1947-03-31 23:00:00 (Mon)
-18000,
0,
'CLT',
    ],
    [
61417627200, #    utc_start 1947-04-01 04:00:00 (Tue)
61422037200, #      utc_end 1947-05-22 05:00:00 (Thu)
61417609200, #  local_start 1947-03-31 23:00:00 (Mon)
61422019200, #    local_end 1947-05-22 00:00:00 (Thu)
-18000,
0,
'CLT',
    ],
    [
61422037200, #    utc_start 1947-05-22 05:00:00 (Thu)
62099064000, #      utc_end 1968-11-03 04:00:00 (Sun)
61422022800, #  local_start 1947-05-22 01:00:00 (Thu)
62099049600, #    local_end 1968-11-03 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62099064000, #    utc_start 1968-11-03 04:00:00 (Sun)
62111761200, #      utc_end 1969-03-30 03:00:00 (Sun)
62099053200, #  local_start 1968-11-03 01:00:00 (Sun)
62111750400, #    local_end 1969-03-30 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62111761200, #    utc_start 1969-03-30 03:00:00 (Sun)
62132328000, #      utc_end 1969-11-23 04:00:00 (Sun)
62111746800, #  local_start 1969-03-29 23:00:00 (Sat)
62132313600, #    local_end 1969-11-23 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62132328000, #    utc_start 1969-11-23 04:00:00 (Sun)
62143210800, #      utc_end 1970-03-29 03:00:00 (Sun)
62132317200, #  local_start 1969-11-23 01:00:00 (Sun)
62143200000, #    local_end 1970-03-29 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62143210800, #    utc_start 1970-03-29 03:00:00 (Sun)
62160148800, #      utc_end 1970-10-11 04:00:00 (Sun)
62143196400, #  local_start 1970-03-28 23:00:00 (Sat)
62160134400, #    local_end 1970-10-11 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62160148800, #    utc_start 1970-10-11 04:00:00 (Sun)
62173450800, #      utc_end 1971-03-14 03:00:00 (Sun)
62160138000, #  local_start 1970-10-11 01:00:00 (Sun)
62173440000, #    local_end 1971-03-14 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62173450800, #    utc_start 1971-03-14 03:00:00 (Sun)
62191598400, #      utc_end 1971-10-10 04:00:00 (Sun)
62173436400, #  local_start 1971-03-13 23:00:00 (Sat)
62191584000, #    local_end 1971-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62191598400, #    utc_start 1971-10-10 04:00:00 (Sun)
62204900400, #      utc_end 1972-03-12 03:00:00 (Sun)
62191587600, #  local_start 1971-10-10 01:00:00 (Sun)
62204889600, #    local_end 1972-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62204900400, #    utc_start 1972-03-12 03:00:00 (Sun)
62223652800, #      utc_end 1972-10-15 04:00:00 (Sun)
62204886000, #  local_start 1972-03-11 23:00:00 (Sat)
62223638400, #    local_end 1972-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62223652800, #    utc_start 1972-10-15 04:00:00 (Sun)
62236350000, #      utc_end 1973-03-11 03:00:00 (Sun)
62223642000, #  local_start 1972-10-15 01:00:00 (Sun)
62236339200, #    local_end 1973-03-11 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62236350000, #    utc_start 1973-03-11 03:00:00 (Sun)
62253892800, #      utc_end 1973-09-30 04:00:00 (Sun)
62236335600, #  local_start 1973-03-10 23:00:00 (Sat)
62253878400, #    local_end 1973-09-30 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62253892800, #    utc_start 1973-09-30 04:00:00 (Sun)
62267799600, #      utc_end 1974-03-10 03:00:00 (Sun)
62253882000, #  local_start 1973-09-30 01:00:00 (Sun)
62267788800, #    local_end 1974-03-10 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62267799600, #    utc_start 1974-03-10 03:00:00 (Sun)
62286552000, #      utc_end 1974-10-13 04:00:00 (Sun)
62267785200, #  local_start 1974-03-09 23:00:00 (Sat)
62286537600, #    local_end 1974-10-13 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62286552000, #    utc_start 1974-10-13 04:00:00 (Sun)
62299249200, #      utc_end 1975-03-09 03:00:00 (Sun)
62286541200, #  local_start 1974-10-13 01:00:00 (Sun)
62299238400, #    local_end 1975-03-09 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62299249200, #    utc_start 1975-03-09 03:00:00 (Sun)
62318001600, #      utc_end 1975-10-12 04:00:00 (Sun)
62299234800, #  local_start 1975-03-08 23:00:00 (Sat)
62317987200, #    local_end 1975-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62318001600, #    utc_start 1975-10-12 04:00:00 (Sun)
62331303600, #      utc_end 1976-03-14 03:00:00 (Sun)
62317990800, #  local_start 1975-10-12 01:00:00 (Sun)
62331292800, #    local_end 1976-03-14 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62331303600, #    utc_start 1976-03-14 03:00:00 (Sun)
62349451200, #      utc_end 1976-10-10 04:00:00 (Sun)
62331289200, #  local_start 1976-03-13 23:00:00 (Sat)
62349436800, #    local_end 1976-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62349451200, #    utc_start 1976-10-10 04:00:00 (Sun)
62362753200, #      utc_end 1977-03-13 03:00:00 (Sun)
62349440400, #  local_start 1976-10-10 01:00:00 (Sun)
62362742400, #    local_end 1977-03-13 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62362753200, #    utc_start 1977-03-13 03:00:00 (Sun)
62380900800, #      utc_end 1977-10-09 04:00:00 (Sun)
62362738800, #  local_start 1977-03-12 23:00:00 (Sat)
62380886400, #    local_end 1977-10-09 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62380900800, #    utc_start 1977-10-09 04:00:00 (Sun)
62394202800, #      utc_end 1978-03-12 03:00:00 (Sun)
62380890000, #  local_start 1977-10-09 01:00:00 (Sun)
62394192000, #    local_end 1978-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62394202800, #    utc_start 1978-03-12 03:00:00 (Sun)
62412955200, #      utc_end 1978-10-15 04:00:00 (Sun)
62394188400, #  local_start 1978-03-11 23:00:00 (Sat)
62412940800, #    local_end 1978-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62412955200, #    utc_start 1978-10-15 04:00:00 (Sun)
62425652400, #      utc_end 1979-03-11 03:00:00 (Sun)
62412944400, #  local_start 1978-10-15 01:00:00 (Sun)
62425641600, #    local_end 1979-03-11 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62425652400, #    utc_start 1979-03-11 03:00:00 (Sun)
62444404800, #      utc_end 1979-10-14 04:00:00 (Sun)
62425638000, #  local_start 1979-03-10 23:00:00 (Sat)
62444390400, #    local_end 1979-10-14 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62444404800, #    utc_start 1979-10-14 04:00:00 (Sun)
62457102000, #      utc_end 1980-03-09 03:00:00 (Sun)
62444394000, #  local_start 1979-10-14 01:00:00 (Sun)
62457091200, #    local_end 1980-03-09 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62457102000, #    utc_start 1980-03-09 03:00:00 (Sun)
62475854400, #      utc_end 1980-10-12 04:00:00 (Sun)
62457087600, #  local_start 1980-03-08 23:00:00 (Sat)
62475840000, #    local_end 1980-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62475854400, #    utc_start 1980-10-12 04:00:00 (Sun)
62489156400, #      utc_end 1981-03-15 03:00:00 (Sun)
62475843600, #  local_start 1980-10-12 01:00:00 (Sun)
62489145600, #    local_end 1981-03-15 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62489156400, #    utc_start 1981-03-15 03:00:00 (Sun)
62507304000, #      utc_end 1981-10-11 04:00:00 (Sun)
62489142000, #  local_start 1981-03-14 23:00:00 (Sat)
62507289600, #    local_end 1981-10-11 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62507304000, #    utc_start 1981-10-11 04:00:00 (Sun)
62520606000, #      utc_end 1982-03-14 03:00:00 (Sun)
62507293200, #  local_start 1981-10-11 01:00:00 (Sun)
62520595200, #    local_end 1982-03-14 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62520606000, #    utc_start 1982-03-14 03:00:00 (Sun)
62538753600, #      utc_end 1982-10-10 04:00:00 (Sun)
62520591600, #  local_start 1982-03-13 23:00:00 (Sat)
62538739200, #    local_end 1982-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62538753600, #    utc_start 1982-10-10 04:00:00 (Sun)
62552055600, #      utc_end 1983-03-13 03:00:00 (Sun)
62538742800, #  local_start 1982-10-10 01:00:00 (Sun)
62552044800, #    local_end 1983-03-13 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62552055600, #    utc_start 1983-03-13 03:00:00 (Sun)
62570203200, #      utc_end 1983-10-09 04:00:00 (Sun)
62552041200, #  local_start 1983-03-12 23:00:00 (Sat)
62570188800, #    local_end 1983-10-09 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62570203200, #    utc_start 1983-10-09 04:00:00 (Sun)
62583505200, #      utc_end 1984-03-11 03:00:00 (Sun)
62570192400, #  local_start 1983-10-09 01:00:00 (Sun)
62583494400, #    local_end 1984-03-11 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62583505200, #    utc_start 1984-03-11 03:00:00 (Sun)
62602257600, #      utc_end 1984-10-14 04:00:00 (Sun)
62583490800, #  local_start 1984-03-10 23:00:00 (Sat)
62602243200, #    local_end 1984-10-14 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62602257600, #    utc_start 1984-10-14 04:00:00 (Sun)
62614954800, #      utc_end 1985-03-10 03:00:00 (Sun)
62602246800, #  local_start 1984-10-14 01:00:00 (Sun)
62614944000, #    local_end 1985-03-10 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62614954800, #    utc_start 1985-03-10 03:00:00 (Sun)
62633707200, #      utc_end 1985-10-13 04:00:00 (Sun)
62614940400, #  local_start 1985-03-09 23:00:00 (Sat)
62633692800, #    local_end 1985-10-13 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62633707200, #    utc_start 1985-10-13 04:00:00 (Sun)
62646404400, #      utc_end 1986-03-09 03:00:00 (Sun)
62633696400, #  local_start 1985-10-13 01:00:00 (Sun)
62646393600, #    local_end 1986-03-09 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62646404400, #    utc_start 1986-03-09 03:00:00 (Sun)
62665156800, #      utc_end 1986-10-12 04:00:00 (Sun)
62646390000, #  local_start 1986-03-08 23:00:00 (Sat)
62665142400, #    local_end 1986-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62665156800, #    utc_start 1986-10-12 04:00:00 (Sun)
62680878000, #      utc_end 1987-04-12 03:00:00 (Sun)
62665146000, #  local_start 1986-10-12 01:00:00 (Sun)
62680867200, #    local_end 1987-04-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62680878000, #    utc_start 1987-04-12 03:00:00 (Sun)
62696606400, #      utc_end 1987-10-11 04:00:00 (Sun)
62680863600, #  local_start 1987-04-11 23:00:00 (Sat)
62696592000, #    local_end 1987-10-11 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62696606400, #    utc_start 1987-10-11 04:00:00 (Sun)
62709908400, #      utc_end 1988-03-13 03:00:00 (Sun)
62696595600, #  local_start 1987-10-11 01:00:00 (Sun)
62709897600, #    local_end 1988-03-13 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62709908400, #    utc_start 1988-03-13 03:00:00 (Sun)
62727451200, #      utc_end 1988-10-02 04:00:00 (Sun)
62709894000, #  local_start 1988-03-12 23:00:00 (Sat)
62727436800, #    local_end 1988-10-02 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62727451200, #    utc_start 1988-10-02 04:00:00 (Sun)
62741358000, #      utc_end 1989-03-12 03:00:00 (Sun)
62727440400, #  local_start 1988-10-02 01:00:00 (Sun)
62741347200, #    local_end 1989-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62741358000, #    utc_start 1989-03-12 03:00:00 (Sun)
62760110400, #      utc_end 1989-10-15 04:00:00 (Sun)
62741343600, #  local_start 1989-03-11 23:00:00 (Sat)
62760096000, #    local_end 1989-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62760110400, #    utc_start 1989-10-15 04:00:00 (Sun)
62773412400, #      utc_end 1990-03-18 03:00:00 (Sun)
62760099600, #  local_start 1989-10-15 01:00:00 (Sun)
62773401600, #    local_end 1990-03-18 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62773412400, #    utc_start 1990-03-18 03:00:00 (Sun)
62789140800, #      utc_end 1990-09-16 04:00:00 (Sun)
62773398000, #  local_start 1990-03-17 23:00:00 (Sat)
62789126400, #    local_end 1990-09-16 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62789140800, #    utc_start 1990-09-16 04:00:00 (Sun)
62804257200, #      utc_end 1991-03-10 03:00:00 (Sun)
62789130000, #  local_start 1990-09-16 01:00:00 (Sun)
62804246400, #    local_end 1991-03-10 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62804257200, #    utc_start 1991-03-10 03:00:00 (Sun)
62823009600, #      utc_end 1991-10-13 04:00:00 (Sun)
62804242800, #  local_start 1991-03-09 23:00:00 (Sat)
62822995200, #    local_end 1991-10-13 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62823009600, #    utc_start 1991-10-13 04:00:00 (Sun)
62836311600, #      utc_end 1992-03-15 03:00:00 (Sun)
62822998800, #  local_start 1991-10-13 01:00:00 (Sun)
62836300800, #    local_end 1992-03-15 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62836311600, #    utc_start 1992-03-15 03:00:00 (Sun)
62854459200, #      utc_end 1992-10-11 04:00:00 (Sun)
62836297200, #  local_start 1992-03-14 23:00:00 (Sat)
62854444800, #    local_end 1992-10-11 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62854459200, #    utc_start 1992-10-11 04:00:00 (Sun)
62867761200, #      utc_end 1993-03-14 03:00:00 (Sun)
62854448400, #  local_start 1992-10-11 01:00:00 (Sun)
62867750400, #    local_end 1993-03-14 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62867761200, #    utc_start 1993-03-14 03:00:00 (Sun)
62885908800, #      utc_end 1993-10-10 04:00:00 (Sun)
62867746800, #  local_start 1993-03-13 23:00:00 (Sat)
62885894400, #    local_end 1993-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62885908800, #    utc_start 1993-10-10 04:00:00 (Sun)
62899210800, #      utc_end 1994-03-13 03:00:00 (Sun)
62885898000, #  local_start 1993-10-10 01:00:00 (Sun)
62899200000, #    local_end 1994-03-13 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62899210800, #    utc_start 1994-03-13 03:00:00 (Sun)
62917358400, #      utc_end 1994-10-09 04:00:00 (Sun)
62899196400, #  local_start 1994-03-12 23:00:00 (Sat)
62917344000, #    local_end 1994-10-09 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62917358400, #    utc_start 1994-10-09 04:00:00 (Sun)
62930660400, #      utc_end 1995-03-12 03:00:00 (Sun)
62917347600, #  local_start 1994-10-09 01:00:00 (Sun)
62930649600, #    local_end 1995-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62930660400, #    utc_start 1995-03-12 03:00:00 (Sun)
62949412800, #      utc_end 1995-10-15 04:00:00 (Sun)
62930646000, #  local_start 1995-03-11 23:00:00 (Sat)
62949398400, #    local_end 1995-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62949412800, #    utc_start 1995-10-15 04:00:00 (Sun)
62962110000, #      utc_end 1996-03-10 03:00:00 (Sun)
62949402000, #  local_start 1995-10-15 01:00:00 (Sun)
62962099200, #    local_end 1996-03-10 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62962110000, #    utc_start 1996-03-10 03:00:00 (Sun)
62980862400, #      utc_end 1996-10-13 04:00:00 (Sun)
62962095600, #  local_start 1996-03-09 23:00:00 (Sat)
62980848000, #    local_end 1996-10-13 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
62980862400, #    utc_start 1996-10-13 04:00:00 (Sun)
62995374000, #      utc_end 1997-03-30 03:00:00 (Sun)
62980851600, #  local_start 1996-10-13 01:00:00 (Sun)
62995363200, #    local_end 1997-03-30 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
62995374000, #    utc_start 1997-03-30 03:00:00 (Sun)
63012312000, #      utc_end 1997-10-12 04:00:00 (Sun)
62995359600, #  local_start 1997-03-29 23:00:00 (Sat)
63012297600, #    local_end 1997-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63012312000, #    utc_start 1997-10-12 04:00:00 (Sun)
63025614000, #      utc_end 1998-03-15 03:00:00 (Sun)
63012301200, #  local_start 1997-10-12 01:00:00 (Sun)
63025603200, #    local_end 1998-03-15 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63025614000, #    utc_start 1998-03-15 03:00:00 (Sun)
63042552000, #      utc_end 1998-09-27 04:00:00 (Sun)
63025599600, #  local_start 1998-03-14 23:00:00 (Sat)
63042537600, #    local_end 1998-09-27 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63042552000, #    utc_start 1998-09-27 04:00:00 (Sun)
63058878000, #      utc_end 1999-04-04 03:00:00 (Sun)
63042541200, #  local_start 1998-09-27 01:00:00 (Sun)
63058867200, #    local_end 1999-04-04 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63058878000, #    utc_start 1999-04-04 03:00:00 (Sun)
63075211200, #      utc_end 1999-10-10 04:00:00 (Sun)
63058863600, #  local_start 1999-04-03 23:00:00 (Sat)
63075196800, #    local_end 1999-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63075211200, #    utc_start 1999-10-10 04:00:00 (Sun)
63088513200, #      utc_end 2000-03-12 03:00:00 (Sun)
63075200400, #  local_start 1999-10-10 01:00:00 (Sun)
63088502400, #    local_end 2000-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63088513200, #    utc_start 2000-03-12 03:00:00 (Sun)
63107265600, #      utc_end 2000-10-15 04:00:00 (Sun)
63088498800, #  local_start 2000-03-11 23:00:00 (Sat)
63107251200, #    local_end 2000-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63107265600, #    utc_start 2000-10-15 04:00:00 (Sun)
63119962800, #      utc_end 2001-03-11 03:00:00 (Sun)
63107254800, #  local_start 2000-10-15 01:00:00 (Sun)
63119952000, #    local_end 2001-03-11 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63119962800, #    utc_start 2001-03-11 03:00:00 (Sun)
63138715200, #      utc_end 2001-10-14 04:00:00 (Sun)
63119948400, #  local_start 2001-03-10 23:00:00 (Sat)
63138700800, #    local_end 2001-10-14 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63138715200, #    utc_start 2001-10-14 04:00:00 (Sun)
63151412400, #      utc_end 2002-03-10 03:00:00 (Sun)
63138704400, #  local_start 2001-10-14 01:00:00 (Sun)
63151401600, #    local_end 2002-03-10 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63151412400, #    utc_start 2002-03-10 03:00:00 (Sun)
63170164800, #      utc_end 2002-10-13 04:00:00 (Sun)
63151398000, #  local_start 2002-03-09 23:00:00 (Sat)
63170150400, #    local_end 2002-10-13 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63170164800, #    utc_start 2002-10-13 04:00:00 (Sun)
63182862000, #      utc_end 2003-03-09 03:00:00 (Sun)
63170154000, #  local_start 2002-10-13 01:00:00 (Sun)
63182851200, #    local_end 2003-03-09 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63182862000, #    utc_start 2003-03-09 03:00:00 (Sun)
63201614400, #      utc_end 2003-10-12 04:00:00 (Sun)
63182847600, #  local_start 2003-03-08 23:00:00 (Sat)
63201600000, #    local_end 2003-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63201614400, #    utc_start 2003-10-12 04:00:00 (Sun)
63214916400, #      utc_end 2004-03-14 03:00:00 (Sun)
63201603600, #  local_start 2003-10-12 01:00:00 (Sun)
63214905600, #    local_end 2004-03-14 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63214916400, #    utc_start 2004-03-14 03:00:00 (Sun)
63233064000, #      utc_end 2004-10-10 04:00:00 (Sun)
63214902000, #  local_start 2004-03-13 23:00:00 (Sat)
63233049600, #    local_end 2004-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63233064000, #    utc_start 2004-10-10 04:00:00 (Sun)
63246366000, #      utc_end 2005-03-13 03:00:00 (Sun)
63233053200, #  local_start 2004-10-10 01:00:00 (Sun)
63246355200, #    local_end 2005-03-13 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63246366000, #    utc_start 2005-03-13 03:00:00 (Sun)
63264513600, #      utc_end 2005-10-09 04:00:00 (Sun)
63246351600, #  local_start 2005-03-12 23:00:00 (Sat)
63264499200, #    local_end 2005-10-09 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63264513600, #    utc_start 2005-10-09 04:00:00 (Sun)
63277815600, #      utc_end 2006-03-12 03:00:00 (Sun)
63264502800, #  local_start 2005-10-09 01:00:00 (Sun)
63277804800, #    local_end 2006-03-12 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63277815600, #    utc_start 2006-03-12 03:00:00 (Sun)
63296568000, #      utc_end 2006-10-15 04:00:00 (Sun)
63277801200, #  local_start 2006-03-11 23:00:00 (Sat)
63296553600, #    local_end 2006-10-15 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63296568000, #    utc_start 2006-10-15 04:00:00 (Sun)
63309265200, #      utc_end 2007-03-11 03:00:00 (Sun)
63296557200, #  local_start 2006-10-15 01:00:00 (Sun)
63309254400, #    local_end 2007-03-11 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63309265200, #    utc_start 2007-03-11 03:00:00 (Sun)
63328017600, #      utc_end 2007-10-14 04:00:00 (Sun)
63309250800, #  local_start 2007-03-10 23:00:00 (Sat)
63328003200, #    local_end 2007-10-14 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63328017600, #    utc_start 2007-10-14 04:00:00 (Sun)
63342529200, #      utc_end 2008-03-30 03:00:00 (Sun)
63328006800, #  local_start 2007-10-14 01:00:00 (Sun)
63342518400, #    local_end 2008-03-30 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63342529200, #    utc_start 2008-03-30 03:00:00 (Sun)
63359467200, #      utc_end 2008-10-12 04:00:00 (Sun)
63342514800, #  local_start 2008-03-29 23:00:00 (Sat)
63359452800, #    local_end 2008-10-12 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63359467200, #    utc_start 2008-10-12 04:00:00 (Sun)
63372769200, #      utc_end 2009-03-15 03:00:00 (Sun)
63359456400, #  local_start 2008-10-12 01:00:00 (Sun)
63372758400, #    local_end 2009-03-15 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63372769200, #    utc_start 2009-03-15 03:00:00 (Sun)
63390916800, #      utc_end 2009-10-11 04:00:00 (Sun)
63372754800, #  local_start 2009-03-14 23:00:00 (Sat)
63390902400, #    local_end 2009-10-11 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63390916800, #    utc_start 2009-10-11 04:00:00 (Sun)
63406033200, #      utc_end 2010-04-04 03:00:00 (Sun)
63390906000, #  local_start 2009-10-11 01:00:00 (Sun)
63406022400, #    local_end 2010-04-04 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63406033200, #    utc_start 2010-04-04 03:00:00 (Sun)
63422366400, #      utc_end 2010-10-10 04:00:00 (Sun)
63406018800, #  local_start 2010-04-03 23:00:00 (Sat)
63422352000, #    local_end 2010-10-10 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63422366400, #    utc_start 2010-10-10 04:00:00 (Sun)
63440506800, #      utc_end 2011-05-08 03:00:00 (Sun)
63422355600, #  local_start 2010-10-10 01:00:00 (Sun)
63440496000, #    local_end 2011-05-08 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63440506800, #    utc_start 2011-05-08 03:00:00 (Sun)
63449582400, #      utc_end 2011-08-21 04:00:00 (Sun)
63440492400, #  local_start 2011-05-07 23:00:00 (Sat)
63449568000, #    local_end 2011-08-21 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63449582400, #    utc_start 2011-08-21 04:00:00 (Sun)
63471351600, #      utc_end 2012-04-29 03:00:00 (Sun)
63449571600, #  local_start 2011-08-21 01:00:00 (Sun)
63471340800, #    local_end 2012-04-29 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63471351600, #    utc_start 2012-04-29 03:00:00 (Sun)
63482241600, #      utc_end 2012-09-02 04:00:00 (Sun)
63471337200, #  local_start 2012-04-28 23:00:00 (Sat)
63482227200, #    local_end 2012-09-02 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63482241600, #    utc_start 2012-09-02 04:00:00 (Sun)
63502801200, #      utc_end 2013-04-28 03:00:00 (Sun)
63482230800, #  local_start 2012-09-02 01:00:00 (Sun)
63502790400, #    local_end 2013-04-28 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63502801200, #    utc_start 2013-04-28 03:00:00 (Sun)
63514296000, #      utc_end 2013-09-08 04:00:00 (Sun)
63502786800, #  local_start 2013-04-27 23:00:00 (Sat)
63514281600, #    local_end 2013-09-08 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63514296000, #    utc_start 2013-09-08 04:00:00 (Sun)
63534250800, #      utc_end 2014-04-27 03:00:00 (Sun)
63514285200, #  local_start 2013-09-08 01:00:00 (Sun)
63534240000, #    local_end 2014-04-27 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63534250800, #    utc_start 2014-04-27 03:00:00 (Sun)
63545745600, #      utc_end 2014-09-07 04:00:00 (Sun)
63534236400, #  local_start 2014-04-26 23:00:00 (Sat)
63545731200, #    local_end 2014-09-07 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63545745600, #    utc_start 2014-09-07 04:00:00 (Sun)
63565700400, #      utc_end 2015-04-26 03:00:00 (Sun)
63545734800, #  local_start 2014-09-07 01:00:00 (Sun)
63565689600, #    local_end 2015-04-26 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63565700400, #    utc_start 2015-04-26 03:00:00 (Sun)
63577195200, #      utc_end 2015-09-06 04:00:00 (Sun)
63565686000, #  local_start 2015-04-25 23:00:00 (Sat)
63577180800, #    local_end 2015-09-06 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63577195200, #    utc_start 2015-09-06 04:00:00 (Sun)
63597150000, #      utc_end 2016-04-24 03:00:00 (Sun)
63577184400, #  local_start 2015-09-06 01:00:00 (Sun)
63597139200, #    local_end 2016-04-24 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63597150000, #    utc_start 2016-04-24 03:00:00 (Sun)
63608644800, #      utc_end 2016-09-04 04:00:00 (Sun)
63597135600, #  local_start 2016-04-23 23:00:00 (Sat)
63608630400, #    local_end 2016-09-04 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63608644800, #    utc_start 2016-09-04 04:00:00 (Sun)
63628599600, #      utc_end 2017-04-23 03:00:00 (Sun)
63608634000, #  local_start 2016-09-04 01:00:00 (Sun)
63628588800, #    local_end 2017-04-23 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63628599600, #    utc_start 2017-04-23 03:00:00 (Sun)
63640094400, #      utc_end 2017-09-03 04:00:00 (Sun)
63628585200, #  local_start 2017-04-22 23:00:00 (Sat)
63640080000, #    local_end 2017-09-03 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63640094400, #    utc_start 2017-09-03 04:00:00 (Sun)
63660654000, #      utc_end 2018-04-29 03:00:00 (Sun)
63640083600, #  local_start 2017-09-03 01:00:00 (Sun)
63660643200, #    local_end 2018-04-29 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63660654000, #    utc_start 2018-04-29 03:00:00 (Sun)
63671544000, #      utc_end 2018-09-02 04:00:00 (Sun)
63660639600, #  local_start 2018-04-28 23:00:00 (Sat)
63671529600, #    local_end 2018-09-02 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63671544000, #    utc_start 2018-09-02 04:00:00 (Sun)
63692103600, #      utc_end 2019-04-28 03:00:00 (Sun)
63671533200, #  local_start 2018-09-02 01:00:00 (Sun)
63692092800, #    local_end 2019-04-28 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63692103600, #    utc_start 2019-04-28 03:00:00 (Sun)
63703598400, #      utc_end 2019-09-08 04:00:00 (Sun)
63692089200, #  local_start 2019-04-27 23:00:00 (Sat)
63703584000, #    local_end 2019-09-08 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63703598400, #    utc_start 2019-09-08 04:00:00 (Sun)
63723553200, #      utc_end 2020-04-26 03:00:00 (Sun)
63703587600, #  local_start 2019-09-08 01:00:00 (Sun)
63723542400, #    local_end 2020-04-26 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63723553200, #    utc_start 2020-04-26 03:00:00 (Sun)
63735048000, #      utc_end 2020-09-06 04:00:00 (Sun)
63723538800, #  local_start 2020-04-25 23:00:00 (Sat)
63735033600, #    local_end 2020-09-06 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63735048000, #    utc_start 2020-09-06 04:00:00 (Sun)
63755002800, #      utc_end 2021-04-25 03:00:00 (Sun)
63735037200, #  local_start 2020-09-06 01:00:00 (Sun)
63754992000, #    local_end 2021-04-25 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63755002800, #    utc_start 2021-04-25 03:00:00 (Sun)
63766497600, #      utc_end 2021-09-05 04:00:00 (Sun)
63754988400, #  local_start 2021-04-24 23:00:00 (Sat)
63766483200, #    local_end 2021-09-05 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63766497600, #    utc_start 2021-09-05 04:00:00 (Sun)
63786452400, #      utc_end 2022-04-24 03:00:00 (Sun)
63766486800, #  local_start 2021-09-05 01:00:00 (Sun)
63786441600, #    local_end 2022-04-24 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63786452400, #    utc_start 2022-04-24 03:00:00 (Sun)
63797947200, #      utc_end 2022-09-04 04:00:00 (Sun)
63786438000, #  local_start 2022-04-23 23:00:00 (Sat)
63797932800, #    local_end 2022-09-04 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63797947200, #    utc_start 2022-09-04 04:00:00 (Sun)
63817902000, #      utc_end 2023-04-23 03:00:00 (Sun)
63797936400, #  local_start 2022-09-04 01:00:00 (Sun)
63817891200, #    local_end 2023-04-23 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63817902000, #    utc_start 2023-04-23 03:00:00 (Sun)
63829396800, #      utc_end 2023-09-03 04:00:00 (Sun)
63817887600, #  local_start 2023-04-22 23:00:00 (Sat)
63829382400, #    local_end 2023-09-03 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
    [
63829396800, #    utc_start 2023-09-03 04:00:00 (Sun)
63849956400, #      utc_end 2024-04-28 03:00:00 (Sun)
63829386000, #  local_start 2023-09-03 01:00:00 (Sun)
63849945600, #    local_end 2024-04-28 00:00:00 (Sun)
-10800,
1,
'CLST',
    ],
    [
63849956400, #    utc_start 2024-04-28 03:00:00 (Sun)
63861451200, #      utc_end 2024-09-08 04:00:00 (Sun)
63849942000, #  local_start 2024-04-27 23:00:00 (Sat)
63861436800, #    local_end 2024-09-08 00:00:00 (Sun)
-14400,
0,
'CLT',
    ],
];

sub olson_version { '2013d' }

sub has_dst_changes { 65 }

sub _max_year { 2023 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}

sub _last_offset { -14400 }

my $last_observance = bless( {
  'format' => 'CL%sT',
  'gmtoff' => '-4:00',
  'local_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 710903,
    'local_rd_secs' => 3600,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 710903,
    'utc_rd_secs' => 3600,
    'utc_year' => 1948
  }, 'DateTime' ),
  'offset_from_std' => 0,
  'offset_from_utc' => -14400,
  'until' => [],
  'utc_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 710903,
    'local_rd_secs' => 18000,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 710903,
    'utc_rd_secs' => 18000,
    'utc_year' => 1948
  }, 'DateTime' )
}, 'DateTime::TimeZone::OlsonDB::Observance' )
;
sub _last_observance { $last_observance }

my $rules = [
  bless( {
    'at' => '4:00u',
    'from' => '2012',
    'in' => 'Sep',
    'letter' => 'S',
    'name' => 'Chile',
    'offset_from_std' => 3600,
    'on' => 'Sun>=2',
    'save' => '1:00',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' ),
  bless( {
    'at' => '3:00u',
    'from' => '2012',
    'in' => 'Apr',
    'letter' => '',
    'name' => 'Chile',
    'offset_from_std' => 0,
    'on' => 'Sun>=23',
    'save' => '0',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' )
]
;
sub _rules { $rules }


1;

