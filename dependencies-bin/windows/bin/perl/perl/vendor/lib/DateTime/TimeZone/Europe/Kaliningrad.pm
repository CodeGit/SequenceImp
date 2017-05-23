# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/6Pwc8w6J1M/europe.  Olson data version 2013d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::Europe::Kaliningrad;
{
  $DateTime::TimeZone::Europe::Kaliningrad::VERSION = '1.60';
}
BEGIN {
  $DateTime::TimeZone::Europe::Kaliningrad::AUTHORITY = 'cpan:DROLSKY';
}

use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::Europe::Kaliningrad::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
59713627080, #      utc_end 1893-03-31 22:38:00 (Fri)
DateTime::TimeZone::NEG_INFINITY, #  local_start
59713632000, #    local_end 1893-04-01 00:00:00 (Sat)
4920,
0,
'LMT',
    ],
    [
59713627080, #    utc_start 1893-03-31 22:38:00 (Fri)
60441976800, #      utc_end 1916-04-30 22:00:00 (Sun)
59713630680, #  local_start 1893-03-31 23:38:00 (Fri)
60441980400, #    local_end 1916-04-30 23:00:00 (Sun)
3600,
0,
'CET',
    ],
    [
60441976800, #    utc_start 1916-04-30 22:00:00 (Sun)
60455199600, #      utc_end 1916-09-30 23:00:00 (Sat)
60441984000, #  local_start 1916-05-01 00:00:00 (Mon)
60455206800, #    local_end 1916-10-01 01:00:00 (Sun)
7200,
1,
'CEST',
    ],
    [
60455199600, #    utc_start 1916-09-30 23:00:00 (Sat)
60472227600, #      utc_end 1917-04-16 01:00:00 (Mon)
60455203200, #  local_start 1916-10-01 00:00:00 (Sun)
60472231200, #    local_end 1917-04-16 02:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
60472227600, #    utc_start 1917-04-16 01:00:00 (Mon)
60485533200, #      utc_end 1917-09-17 01:00:00 (Mon)
60472234800, #  local_start 1917-04-16 03:00:00 (Mon)
60485540400, #    local_end 1917-09-17 03:00:00 (Mon)
7200,
1,
'CEST',
    ],
    [
60485533200, #    utc_start 1917-09-17 01:00:00 (Mon)
60503677200, #      utc_end 1918-04-15 01:00:00 (Mon)
60485536800, #  local_start 1917-09-17 02:00:00 (Mon)
60503680800, #    local_end 1918-04-15 02:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
60503677200, #    utc_start 1918-04-15 01:00:00 (Mon)
60516982800, #      utc_end 1918-09-16 01:00:00 (Mon)
60503684400, #  local_start 1918-04-15 03:00:00 (Mon)
60516990000, #    local_end 1918-09-16 03:00:00 (Mon)
7200,
1,
'CEST',
    ],
    [
60516982800, #    utc_start 1918-09-16 01:00:00 (Mon)
61196778000, #      utc_end 1940-04-01 01:00:00 (Mon)
60516986400, #  local_start 1918-09-16 02:00:00 (Mon)
61196781600, #    local_end 1940-04-01 02:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
61196778000, #    utc_start 1940-04-01 01:00:00 (Mon)
61278426000, #      utc_end 1942-11-02 01:00:00 (Mon)
61196785200, #  local_start 1940-04-01 03:00:00 (Mon)
61278433200, #    local_end 1942-11-02 03:00:00 (Mon)
7200,
1,
'CEST',
    ],
    [
61278426000, #    utc_start 1942-11-02 01:00:00 (Mon)
61291126800, #      utc_end 1943-03-29 01:00:00 (Mon)
61278429600, #  local_start 1942-11-02 02:00:00 (Mon)
61291130400, #    local_end 1943-03-29 02:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
61291126800, #    utc_start 1943-03-29 01:00:00 (Mon)
61307456400, #      utc_end 1943-10-04 01:00:00 (Mon)
61291134000, #  local_start 1943-03-29 03:00:00 (Mon)
61307463600, #    local_end 1943-10-04 03:00:00 (Mon)
7200,
1,
'CEST',
    ],
    [
61307456400, #    utc_start 1943-10-04 01:00:00 (Mon)
61323181200, #      utc_end 1944-04-03 01:00:00 (Mon)
61307460000, #  local_start 1943-10-04 02:00:00 (Mon)
61323184800, #    local_end 1944-04-03 02:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
61323181200, #    utc_start 1944-04-03 01:00:00 (Mon)
61338906000, #      utc_end 1944-10-02 01:00:00 (Mon)
61323188400, #  local_start 1944-04-03 03:00:00 (Mon)
61338913200, #    local_end 1944-10-02 03:00:00 (Mon)
7200,
1,
'CEST',
    ],
    [
61338906000, #    utc_start 1944-10-02 01:00:00 (Mon)
61346761200, #      utc_end 1944-12-31 23:00:00 (Sun)
61338909600, #  local_start 1944-10-02 02:00:00 (Mon)
61346764800, #    local_end 1945-01-01 00:00:00 (Mon)
3600,
0,
'CET',
    ],
    [
61346761200, #    utc_start 1944-12-31 23:00:00 (Sun)
61356952800, #      utc_end 1945-04-28 22:00:00 (Sat)
61346768400, #  local_start 1945-01-01 01:00:00 (Mon)
61356960000, #    local_end 1945-04-29 00:00:00 (Sun)
7200,
0,
'CET',
    ],
    [
61356952800, #    utc_start 1945-04-28 22:00:00 (Sat)
61373019600, #      utc_end 1945-10-31 21:00:00 (Wed)
61356963600, #  local_start 1945-04-29 01:00:00 (Sun)
61373030400, #    local_end 1945-11-01 00:00:00 (Thu)
10800,
1,
'CEST',
    ],
    [
61373019600, #    utc_start 1945-10-31 21:00:00 (Wed)
61378293600, #      utc_end 1945-12-31 22:00:00 (Mon)
61373026800, #  local_start 1945-10-31 23:00:00 (Wed)
61378300800, #    local_end 1946-01-01 00:00:00 (Tue)
7200,
0,
'CET',
    ],
    [
61378293600, #    utc_start 1945-12-31 22:00:00 (Mon)
62490603600, #      utc_end 1981-03-31 21:00:00 (Tue)
61378304400, #  local_start 1946-01-01 01:00:00 (Tue)
62490614400, #    local_end 1981-04-01 00:00:00 (Wed)
10800,
0,
'MSK',
    ],
    [
62490603600, #    utc_start 1981-03-31 21:00:00 (Tue)
62506411200, #      utc_end 1981-09-30 20:00:00 (Wed)
62490618000, #  local_start 1981-04-01 01:00:00 (Wed)
62506425600, #    local_end 1981-10-01 00:00:00 (Thu)
14400,
1,
'MSD',
    ],
    [
62506411200, #    utc_start 1981-09-30 20:00:00 (Wed)
62522139600, #      utc_end 1982-03-31 21:00:00 (Wed)
62506422000, #  local_start 1981-09-30 23:00:00 (Wed)
62522150400, #    local_end 1982-04-01 00:00:00 (Thu)
10800,
0,
'MSK',
    ],
    [
62522139600, #    utc_start 1982-03-31 21:00:00 (Wed)
62537947200, #      utc_end 1982-09-30 20:00:00 (Thu)
62522154000, #  local_start 1982-04-01 01:00:00 (Thu)
62537961600, #    local_end 1982-10-01 00:00:00 (Fri)
14400,
1,
'MSD',
    ],
    [
62537947200, #    utc_start 1982-09-30 20:00:00 (Thu)
62553675600, #      utc_end 1983-03-31 21:00:00 (Thu)
62537958000, #  local_start 1982-09-30 23:00:00 (Thu)
62553686400, #    local_end 1983-04-01 00:00:00 (Fri)
10800,
0,
'MSK',
    ],
    [
62553675600, #    utc_start 1983-03-31 21:00:00 (Thu)
62569483200, #      utc_end 1983-09-30 20:00:00 (Fri)
62553690000, #  local_start 1983-04-01 01:00:00 (Fri)
62569497600, #    local_end 1983-10-01 00:00:00 (Sat)
14400,
1,
'MSD',
    ],
    [
62569483200, #    utc_start 1983-09-30 20:00:00 (Fri)
62585298000, #      utc_end 1984-03-31 21:00:00 (Sat)
62569494000, #  local_start 1983-09-30 23:00:00 (Fri)
62585308800, #    local_end 1984-04-01 00:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62585298000, #    utc_start 1984-03-31 21:00:00 (Sat)
62601030000, #      utc_end 1984-09-29 23:00:00 (Sat)
62585312400, #  local_start 1984-04-01 01:00:00 (Sun)
62601044400, #    local_end 1984-09-30 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62601030000, #    utc_start 1984-09-29 23:00:00 (Sat)
62616754800, #      utc_end 1985-03-30 23:00:00 (Sat)
62601040800, #  local_start 1984-09-30 02:00:00 (Sun)
62616765600, #    local_end 1985-03-31 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62616754800, #    utc_start 1985-03-30 23:00:00 (Sat)
62632479600, #      utc_end 1985-09-28 23:00:00 (Sat)
62616769200, #  local_start 1985-03-31 03:00:00 (Sun)
62632494000, #    local_end 1985-09-29 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62632479600, #    utc_start 1985-09-28 23:00:00 (Sat)
62648204400, #      utc_end 1986-03-29 23:00:00 (Sat)
62632490400, #  local_start 1985-09-29 02:00:00 (Sun)
62648215200, #    local_end 1986-03-30 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62648204400, #    utc_start 1986-03-29 23:00:00 (Sat)
62663929200, #      utc_end 1986-09-27 23:00:00 (Sat)
62648218800, #  local_start 1986-03-30 03:00:00 (Sun)
62663943600, #    local_end 1986-09-28 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62663929200, #    utc_start 1986-09-27 23:00:00 (Sat)
62679654000, #      utc_end 1987-03-28 23:00:00 (Sat)
62663940000, #  local_start 1986-09-28 02:00:00 (Sun)
62679664800, #    local_end 1987-03-29 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62679654000, #    utc_start 1987-03-28 23:00:00 (Sat)
62695378800, #      utc_end 1987-09-26 23:00:00 (Sat)
62679668400, #  local_start 1987-03-29 03:00:00 (Sun)
62695393200, #    local_end 1987-09-27 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62695378800, #    utc_start 1987-09-26 23:00:00 (Sat)
62711103600, #      utc_end 1988-03-26 23:00:00 (Sat)
62695389600, #  local_start 1987-09-27 02:00:00 (Sun)
62711114400, #    local_end 1988-03-27 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62711103600, #    utc_start 1988-03-26 23:00:00 (Sat)
62726828400, #      utc_end 1988-09-24 23:00:00 (Sat)
62711118000, #  local_start 1988-03-27 03:00:00 (Sun)
62726842800, #    local_end 1988-09-25 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62726828400, #    utc_start 1988-09-24 23:00:00 (Sat)
62742553200, #      utc_end 1989-03-25 23:00:00 (Sat)
62726839200, #  local_start 1988-09-25 02:00:00 (Sun)
62742564000, #    local_end 1989-03-26 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62742553200, #    utc_start 1989-03-25 23:00:00 (Sat)
62758278000, #      utc_end 1989-09-23 23:00:00 (Sat)
62742567600, #  local_start 1989-03-26 03:00:00 (Sun)
62758292400, #    local_end 1989-09-24 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62758278000, #    utc_start 1989-09-23 23:00:00 (Sat)
62774002800, #      utc_end 1990-03-24 23:00:00 (Sat)
62758288800, #  local_start 1989-09-24 02:00:00 (Sun)
62774013600, #    local_end 1990-03-25 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62774002800, #    utc_start 1990-03-24 23:00:00 (Sat)
62790332400, #      utc_end 1990-09-29 23:00:00 (Sat)
62774017200, #  local_start 1990-03-25 03:00:00 (Sun)
62790346800, #    local_end 1990-09-30 03:00:00 (Sun)
14400,
1,
'MSD',
    ],
    [
62790332400, #    utc_start 1990-09-29 23:00:00 (Sat)
62806057200, #      utc_end 1991-03-30 23:00:00 (Sat)
62790343200, #  local_start 1990-09-30 02:00:00 (Sun)
62806068000, #    local_end 1991-03-31 02:00:00 (Sun)
10800,
0,
'MSK',
    ],
    [
62806057200, #    utc_start 1991-03-30 23:00:00 (Sat)
62821785600, #      utc_end 1991-09-29 00:00:00 (Sun)
62806068000, #  local_start 1991-03-31 02:00:00 (Sun)
62821796400, #    local_end 1991-09-29 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
62821785600, #    utc_start 1991-09-29 00:00:00 (Sun)
62837499600, #      utc_end 1992-03-28 21:00:00 (Sat)
62821792800, #  local_start 1991-09-29 02:00:00 (Sun)
62837506800, #    local_end 1992-03-28 23:00:00 (Sat)
7200,
0,
'EET',
    ],
    [
62837499600, #    utc_start 1992-03-28 21:00:00 (Sat)
62853220800, #      utc_end 1992-09-26 20:00:00 (Sat)
62837510400, #  local_start 1992-03-29 00:00:00 (Sun)
62853231600, #    local_end 1992-09-26 23:00:00 (Sat)
10800,
1,
'EEST',
    ],
    [
62853220800, #    utc_start 1992-09-26 20:00:00 (Sat)
62868960000, #      utc_end 1993-03-28 00:00:00 (Sun)
62853228000, #  local_start 1992-09-26 22:00:00 (Sat)
62868967200, #    local_end 1993-03-28 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
62868960000, #    utc_start 1993-03-28 00:00:00 (Sun)
62884684800, #      utc_end 1993-09-26 00:00:00 (Sun)
62868970800, #  local_start 1993-03-28 03:00:00 (Sun)
62884695600, #    local_end 1993-09-26 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
62884684800, #    utc_start 1993-09-26 00:00:00 (Sun)
62900409600, #      utc_end 1994-03-27 00:00:00 (Sun)
62884692000, #  local_start 1993-09-26 02:00:00 (Sun)
62900416800, #    local_end 1994-03-27 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
62900409600, #    utc_start 1994-03-27 00:00:00 (Sun)
62916134400, #      utc_end 1994-09-25 00:00:00 (Sun)
62900420400, #  local_start 1994-03-27 03:00:00 (Sun)
62916145200, #    local_end 1994-09-25 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
62916134400, #    utc_start 1994-09-25 00:00:00 (Sun)
62931859200, #      utc_end 1995-03-26 00:00:00 (Sun)
62916141600, #  local_start 1994-09-25 02:00:00 (Sun)
62931866400, #    local_end 1995-03-26 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
62931859200, #    utc_start 1995-03-26 00:00:00 (Sun)
62947584000, #      utc_end 1995-09-24 00:00:00 (Sun)
62931870000, #  local_start 1995-03-26 03:00:00 (Sun)
62947594800, #    local_end 1995-09-24 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
62947584000, #    utc_start 1995-09-24 00:00:00 (Sun)
62963913600, #      utc_end 1996-03-31 00:00:00 (Sun)
62947591200, #  local_start 1995-09-24 02:00:00 (Sun)
62963920800, #    local_end 1996-03-31 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
62963913600, #    utc_start 1996-03-31 00:00:00 (Sun)
62982057600, #      utc_end 1996-10-27 00:00:00 (Sun)
62963924400, #  local_start 1996-03-31 03:00:00 (Sun)
62982068400, #    local_end 1996-10-27 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
62982057600, #    utc_start 1996-10-27 00:00:00 (Sun)
62995363200, #      utc_end 1997-03-30 00:00:00 (Sun)
62982064800, #  local_start 1996-10-27 02:00:00 (Sun)
62995370400, #    local_end 1997-03-30 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
62995363200, #    utc_start 1997-03-30 00:00:00 (Sun)
63013507200, #      utc_end 1997-10-26 00:00:00 (Sun)
62995374000, #  local_start 1997-03-30 03:00:00 (Sun)
63013518000, #    local_end 1997-10-26 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63013507200, #    utc_start 1997-10-26 00:00:00 (Sun)
63026812800, #      utc_end 1998-03-29 00:00:00 (Sun)
63013514400, #  local_start 1997-10-26 02:00:00 (Sun)
63026820000, #    local_end 1998-03-29 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63026812800, #    utc_start 1998-03-29 00:00:00 (Sun)
63044956800, #      utc_end 1998-10-25 00:00:00 (Sun)
63026823600, #  local_start 1998-03-29 03:00:00 (Sun)
63044967600, #    local_end 1998-10-25 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63044956800, #    utc_start 1998-10-25 00:00:00 (Sun)
63058262400, #      utc_end 1999-03-28 00:00:00 (Sun)
63044964000, #  local_start 1998-10-25 02:00:00 (Sun)
63058269600, #    local_end 1999-03-28 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63058262400, #    utc_start 1999-03-28 00:00:00 (Sun)
63077011200, #      utc_end 1999-10-31 00:00:00 (Sun)
63058273200, #  local_start 1999-03-28 03:00:00 (Sun)
63077022000, #    local_end 1999-10-31 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63077011200, #    utc_start 1999-10-31 00:00:00 (Sun)
63089712000, #      utc_end 2000-03-26 00:00:00 (Sun)
63077018400, #  local_start 1999-10-31 02:00:00 (Sun)
63089719200, #    local_end 2000-03-26 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63089712000, #    utc_start 2000-03-26 00:00:00 (Sun)
63108460800, #      utc_end 2000-10-29 00:00:00 (Sun)
63089722800, #  local_start 2000-03-26 03:00:00 (Sun)
63108471600, #    local_end 2000-10-29 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63108460800, #    utc_start 2000-10-29 00:00:00 (Sun)
63121161600, #      utc_end 2001-03-25 00:00:00 (Sun)
63108468000, #  local_start 2000-10-29 02:00:00 (Sun)
63121168800, #    local_end 2001-03-25 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63121161600, #    utc_start 2001-03-25 00:00:00 (Sun)
63139910400, #      utc_end 2001-10-28 00:00:00 (Sun)
63121172400, #  local_start 2001-03-25 03:00:00 (Sun)
63139921200, #    local_end 2001-10-28 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63139910400, #    utc_start 2001-10-28 00:00:00 (Sun)
63153216000, #      utc_end 2002-03-31 00:00:00 (Sun)
63139917600, #  local_start 2001-10-28 02:00:00 (Sun)
63153223200, #    local_end 2002-03-31 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63153216000, #    utc_start 2002-03-31 00:00:00 (Sun)
63171360000, #      utc_end 2002-10-27 00:00:00 (Sun)
63153226800, #  local_start 2002-03-31 03:00:00 (Sun)
63171370800, #    local_end 2002-10-27 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63171360000, #    utc_start 2002-10-27 00:00:00 (Sun)
63184665600, #      utc_end 2003-03-30 00:00:00 (Sun)
63171367200, #  local_start 2002-10-27 02:00:00 (Sun)
63184672800, #    local_end 2003-03-30 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63184665600, #    utc_start 2003-03-30 00:00:00 (Sun)
63202809600, #      utc_end 2003-10-26 00:00:00 (Sun)
63184676400, #  local_start 2003-03-30 03:00:00 (Sun)
63202820400, #    local_end 2003-10-26 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63202809600, #    utc_start 2003-10-26 00:00:00 (Sun)
63216115200, #      utc_end 2004-03-28 00:00:00 (Sun)
63202816800, #  local_start 2003-10-26 02:00:00 (Sun)
63216122400, #    local_end 2004-03-28 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63216115200, #    utc_start 2004-03-28 00:00:00 (Sun)
63234864000, #      utc_end 2004-10-31 00:00:00 (Sun)
63216126000, #  local_start 2004-03-28 03:00:00 (Sun)
63234874800, #    local_end 2004-10-31 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63234864000, #    utc_start 2004-10-31 00:00:00 (Sun)
63247564800, #      utc_end 2005-03-27 00:00:00 (Sun)
63234871200, #  local_start 2004-10-31 02:00:00 (Sun)
63247572000, #    local_end 2005-03-27 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63247564800, #    utc_start 2005-03-27 00:00:00 (Sun)
63266313600, #      utc_end 2005-10-30 00:00:00 (Sun)
63247575600, #  local_start 2005-03-27 03:00:00 (Sun)
63266324400, #    local_end 2005-10-30 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63266313600, #    utc_start 2005-10-30 00:00:00 (Sun)
63279014400, #      utc_end 2006-03-26 00:00:00 (Sun)
63266320800, #  local_start 2005-10-30 02:00:00 (Sun)
63279021600, #    local_end 2006-03-26 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63279014400, #    utc_start 2006-03-26 00:00:00 (Sun)
63297763200, #      utc_end 2006-10-29 00:00:00 (Sun)
63279025200, #  local_start 2006-03-26 03:00:00 (Sun)
63297774000, #    local_end 2006-10-29 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63297763200, #    utc_start 2006-10-29 00:00:00 (Sun)
63310464000, #      utc_end 2007-03-25 00:00:00 (Sun)
63297770400, #  local_start 2006-10-29 02:00:00 (Sun)
63310471200, #    local_end 2007-03-25 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63310464000, #    utc_start 2007-03-25 00:00:00 (Sun)
63329212800, #      utc_end 2007-10-28 00:00:00 (Sun)
63310474800, #  local_start 2007-03-25 03:00:00 (Sun)
63329223600, #    local_end 2007-10-28 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63329212800, #    utc_start 2007-10-28 00:00:00 (Sun)
63342518400, #      utc_end 2008-03-30 00:00:00 (Sun)
63329220000, #  local_start 2007-10-28 02:00:00 (Sun)
63342525600, #    local_end 2008-03-30 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63342518400, #    utc_start 2008-03-30 00:00:00 (Sun)
63360662400, #      utc_end 2008-10-26 00:00:00 (Sun)
63342529200, #  local_start 2008-03-30 03:00:00 (Sun)
63360673200, #    local_end 2008-10-26 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63360662400, #    utc_start 2008-10-26 00:00:00 (Sun)
63373968000, #      utc_end 2009-03-29 00:00:00 (Sun)
63360669600, #  local_start 2008-10-26 02:00:00 (Sun)
63373975200, #    local_end 2009-03-29 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63373968000, #    utc_start 2009-03-29 00:00:00 (Sun)
63392112000, #      utc_end 2009-10-25 00:00:00 (Sun)
63373978800, #  local_start 2009-03-29 03:00:00 (Sun)
63392122800, #    local_end 2009-10-25 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63392112000, #    utc_start 2009-10-25 00:00:00 (Sun)
63405417600, #      utc_end 2010-03-28 00:00:00 (Sun)
63392119200, #  local_start 2009-10-25 02:00:00 (Sun)
63405424800, #    local_end 2010-03-28 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63405417600, #    utc_start 2010-03-28 00:00:00 (Sun)
63424166400, #      utc_end 2010-10-31 00:00:00 (Sun)
63405428400, #  local_start 2010-03-28 03:00:00 (Sun)
63424177200, #    local_end 2010-10-31 03:00:00 (Sun)
10800,
1,
'EEST',
    ],
    [
63424166400, #    utc_start 2010-10-31 00:00:00 (Sun)
63436867200, #      utc_end 2011-03-27 00:00:00 (Sun)
63424173600, #  local_start 2010-10-31 02:00:00 (Sun)
63436874400, #    local_end 2011-03-27 02:00:00 (Sun)
7200,
0,
'EET',
    ],
    [
63436867200, #    utc_start 2011-03-27 00:00:00 (Sun)
DateTime::TimeZone::INFINITY, #      utc_end
63436878000, #  local_start 2011-03-27 03:00:00 (Sun)
DateTime::TimeZone::INFINITY, #    local_end
10800,
0,
'FET',
    ],
];

sub olson_version { '2013d' }

sub has_dst_changes { 37 }

sub _max_year { 2023 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}



1;

