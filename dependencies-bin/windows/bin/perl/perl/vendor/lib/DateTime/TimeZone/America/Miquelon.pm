# This file is auto-generated by the Perl DateTime Suite time zone
# code generator (0.07) This code generator comes with the
# DateTime::TimeZone module distribution in the tools/ directory

#
# Generated from /tmp/6Pwc8w6J1M/northamerica.  Olson data version 2013d
#
# Do not edit this file directly.
#
package DateTime::TimeZone::America::Miquelon;
{
  $DateTime::TimeZone::America::Miquelon::VERSION = '1.60';
}
BEGIN {
  $DateTime::TimeZone::America::Miquelon::AUTHORITY = 'cpan:DROLSKY';
}

use strict;

use Class::Singleton 1.03;
use DateTime::TimeZone;
use DateTime::TimeZone::OlsonDB;

@DateTime::TimeZone::America::Miquelon::ISA = ( 'Class::Singleton', 'DateTime::TimeZone' );

my $spans =
[
    [
DateTime::TimeZone::NEG_INFINITY, #    utc_start
60285354280, #      utc_end 1911-05-15 03:44:40 (Mon)
DateTime::TimeZone::NEG_INFINITY, #  local_start
60285340800, #    local_end 1911-05-15 00:00:00 (Mon)
-13480,
0,
'LMT',
    ],
    [
60285354280, #    utc_start 1911-05-15 03:44:40 (Mon)
62461684800, #      utc_end 1980-05-01 04:00:00 (Thu)
60285339880, #  local_start 1911-05-14 23:44:40 (Sun)
62461670400, #    local_end 1980-05-01 00:00:00 (Thu)
-14400,
0,
'AST',
    ],
    [
62461684800, #    utc_start 1980-05-01 04:00:00 (Thu)
62672151600, #      utc_end 1987-01-01 03:00:00 (Thu)
62461674000, #  local_start 1980-05-01 01:00:00 (Thu)
62672140800, #    local_end 1987-01-01 00:00:00 (Thu)
-10800,
0,
'PMST',
    ],
    [
62672151600, #    utc_start 1987-01-01 03:00:00 (Thu)
62680280400, #      utc_end 1987-04-05 05:00:00 (Sun)
62672140800, #  local_start 1987-01-01 00:00:00 (Thu)
62680269600, #    local_end 1987-04-05 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62680280400, #    utc_start 1987-04-05 05:00:00 (Sun)
62697816000, #      utc_end 1987-10-25 04:00:00 (Sun)
62680273200, #  local_start 1987-04-05 03:00:00 (Sun)
62697808800, #    local_end 1987-10-25 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62697816000, #    utc_start 1987-10-25 04:00:00 (Sun)
62711730000, #      utc_end 1988-04-03 05:00:00 (Sun)
62697805200, #  local_start 1987-10-25 01:00:00 (Sun)
62711719200, #    local_end 1988-04-03 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62711730000, #    utc_start 1988-04-03 05:00:00 (Sun)
62729870400, #      utc_end 1988-10-30 04:00:00 (Sun)
62711722800, #  local_start 1988-04-03 03:00:00 (Sun)
62729863200, #    local_end 1988-10-30 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62729870400, #    utc_start 1988-10-30 04:00:00 (Sun)
62743179600, #      utc_end 1989-04-02 05:00:00 (Sun)
62729859600, #  local_start 1988-10-30 01:00:00 (Sun)
62743168800, #    local_end 1989-04-02 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62743179600, #    utc_start 1989-04-02 05:00:00 (Sun)
62761320000, #      utc_end 1989-10-29 04:00:00 (Sun)
62743172400, #  local_start 1989-04-02 03:00:00 (Sun)
62761312800, #    local_end 1989-10-29 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62761320000, #    utc_start 1989-10-29 04:00:00 (Sun)
62774629200, #      utc_end 1990-04-01 05:00:00 (Sun)
62761309200, #  local_start 1989-10-29 01:00:00 (Sun)
62774618400, #    local_end 1990-04-01 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62774629200, #    utc_start 1990-04-01 05:00:00 (Sun)
62792769600, #      utc_end 1990-10-28 04:00:00 (Sun)
62774622000, #  local_start 1990-04-01 03:00:00 (Sun)
62792762400, #    local_end 1990-10-28 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62792769600, #    utc_start 1990-10-28 04:00:00 (Sun)
62806683600, #      utc_end 1991-04-07 05:00:00 (Sun)
62792758800, #  local_start 1990-10-28 01:00:00 (Sun)
62806672800, #    local_end 1991-04-07 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62806683600, #    utc_start 1991-04-07 05:00:00 (Sun)
62824219200, #      utc_end 1991-10-27 04:00:00 (Sun)
62806676400, #  local_start 1991-04-07 03:00:00 (Sun)
62824212000, #    local_end 1991-10-27 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62824219200, #    utc_start 1991-10-27 04:00:00 (Sun)
62838133200, #      utc_end 1992-04-05 05:00:00 (Sun)
62824208400, #  local_start 1991-10-27 01:00:00 (Sun)
62838122400, #    local_end 1992-04-05 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62838133200, #    utc_start 1992-04-05 05:00:00 (Sun)
62855668800, #      utc_end 1992-10-25 04:00:00 (Sun)
62838126000, #  local_start 1992-04-05 03:00:00 (Sun)
62855661600, #    local_end 1992-10-25 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62855668800, #    utc_start 1992-10-25 04:00:00 (Sun)
62869582800, #      utc_end 1993-04-04 05:00:00 (Sun)
62855658000, #  local_start 1992-10-25 01:00:00 (Sun)
62869572000, #    local_end 1993-04-04 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62869582800, #    utc_start 1993-04-04 05:00:00 (Sun)
62887723200, #      utc_end 1993-10-31 04:00:00 (Sun)
62869575600, #  local_start 1993-04-04 03:00:00 (Sun)
62887716000, #    local_end 1993-10-31 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62887723200, #    utc_start 1993-10-31 04:00:00 (Sun)
62901032400, #      utc_end 1994-04-03 05:00:00 (Sun)
62887712400, #  local_start 1993-10-31 01:00:00 (Sun)
62901021600, #    local_end 1994-04-03 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62901032400, #    utc_start 1994-04-03 05:00:00 (Sun)
62919172800, #      utc_end 1994-10-30 04:00:00 (Sun)
62901025200, #  local_start 1994-04-03 03:00:00 (Sun)
62919165600, #    local_end 1994-10-30 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62919172800, #    utc_start 1994-10-30 04:00:00 (Sun)
62932482000, #      utc_end 1995-04-02 05:00:00 (Sun)
62919162000, #  local_start 1994-10-30 01:00:00 (Sun)
62932471200, #    local_end 1995-04-02 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62932482000, #    utc_start 1995-04-02 05:00:00 (Sun)
62950622400, #      utc_end 1995-10-29 04:00:00 (Sun)
62932474800, #  local_start 1995-04-02 03:00:00 (Sun)
62950615200, #    local_end 1995-10-29 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62950622400, #    utc_start 1995-10-29 04:00:00 (Sun)
62964536400, #      utc_end 1996-04-07 05:00:00 (Sun)
62950611600, #  local_start 1995-10-29 01:00:00 (Sun)
62964525600, #    local_end 1996-04-07 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62964536400, #    utc_start 1996-04-07 05:00:00 (Sun)
62982072000, #      utc_end 1996-10-27 04:00:00 (Sun)
62964529200, #  local_start 1996-04-07 03:00:00 (Sun)
62982064800, #    local_end 1996-10-27 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
62982072000, #    utc_start 1996-10-27 04:00:00 (Sun)
62995986000, #      utc_end 1997-04-06 05:00:00 (Sun)
62982061200, #  local_start 1996-10-27 01:00:00 (Sun)
62995975200, #    local_end 1997-04-06 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
62995986000, #    utc_start 1997-04-06 05:00:00 (Sun)
63013521600, #      utc_end 1997-10-26 04:00:00 (Sun)
62995978800, #  local_start 1997-04-06 03:00:00 (Sun)
63013514400, #    local_end 1997-10-26 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63013521600, #    utc_start 1997-10-26 04:00:00 (Sun)
63027435600, #      utc_end 1998-04-05 05:00:00 (Sun)
63013510800, #  local_start 1997-10-26 01:00:00 (Sun)
63027424800, #    local_end 1998-04-05 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63027435600, #    utc_start 1998-04-05 05:00:00 (Sun)
63044971200, #      utc_end 1998-10-25 04:00:00 (Sun)
63027428400, #  local_start 1998-04-05 03:00:00 (Sun)
63044964000, #    local_end 1998-10-25 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63044971200, #    utc_start 1998-10-25 04:00:00 (Sun)
63058885200, #      utc_end 1999-04-04 05:00:00 (Sun)
63044960400, #  local_start 1998-10-25 01:00:00 (Sun)
63058874400, #    local_end 1999-04-04 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63058885200, #    utc_start 1999-04-04 05:00:00 (Sun)
63077025600, #      utc_end 1999-10-31 04:00:00 (Sun)
63058878000, #  local_start 1999-04-04 03:00:00 (Sun)
63077018400, #    local_end 1999-10-31 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63077025600, #    utc_start 1999-10-31 04:00:00 (Sun)
63090334800, #      utc_end 2000-04-02 05:00:00 (Sun)
63077014800, #  local_start 1999-10-31 01:00:00 (Sun)
63090324000, #    local_end 2000-04-02 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63090334800, #    utc_start 2000-04-02 05:00:00 (Sun)
63108475200, #      utc_end 2000-10-29 04:00:00 (Sun)
63090327600, #  local_start 2000-04-02 03:00:00 (Sun)
63108468000, #    local_end 2000-10-29 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63108475200, #    utc_start 2000-10-29 04:00:00 (Sun)
63121784400, #      utc_end 2001-04-01 05:00:00 (Sun)
63108464400, #  local_start 2000-10-29 01:00:00 (Sun)
63121773600, #    local_end 2001-04-01 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63121784400, #    utc_start 2001-04-01 05:00:00 (Sun)
63139924800, #      utc_end 2001-10-28 04:00:00 (Sun)
63121777200, #  local_start 2001-04-01 03:00:00 (Sun)
63139917600, #    local_end 2001-10-28 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63139924800, #    utc_start 2001-10-28 04:00:00 (Sun)
63153838800, #      utc_end 2002-04-07 05:00:00 (Sun)
63139914000, #  local_start 2001-10-28 01:00:00 (Sun)
63153828000, #    local_end 2002-04-07 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63153838800, #    utc_start 2002-04-07 05:00:00 (Sun)
63171374400, #      utc_end 2002-10-27 04:00:00 (Sun)
63153831600, #  local_start 2002-04-07 03:00:00 (Sun)
63171367200, #    local_end 2002-10-27 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63171374400, #    utc_start 2002-10-27 04:00:00 (Sun)
63185288400, #      utc_end 2003-04-06 05:00:00 (Sun)
63171363600, #  local_start 2002-10-27 01:00:00 (Sun)
63185277600, #    local_end 2003-04-06 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63185288400, #    utc_start 2003-04-06 05:00:00 (Sun)
63202824000, #      utc_end 2003-10-26 04:00:00 (Sun)
63185281200, #  local_start 2003-04-06 03:00:00 (Sun)
63202816800, #    local_end 2003-10-26 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63202824000, #    utc_start 2003-10-26 04:00:00 (Sun)
63216738000, #      utc_end 2004-04-04 05:00:00 (Sun)
63202813200, #  local_start 2003-10-26 01:00:00 (Sun)
63216727200, #    local_end 2004-04-04 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63216738000, #    utc_start 2004-04-04 05:00:00 (Sun)
63234878400, #      utc_end 2004-10-31 04:00:00 (Sun)
63216730800, #  local_start 2004-04-04 03:00:00 (Sun)
63234871200, #    local_end 2004-10-31 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63234878400, #    utc_start 2004-10-31 04:00:00 (Sun)
63248187600, #      utc_end 2005-04-03 05:00:00 (Sun)
63234867600, #  local_start 2004-10-31 01:00:00 (Sun)
63248176800, #    local_end 2005-04-03 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63248187600, #    utc_start 2005-04-03 05:00:00 (Sun)
63266328000, #      utc_end 2005-10-30 04:00:00 (Sun)
63248180400, #  local_start 2005-04-03 03:00:00 (Sun)
63266320800, #    local_end 2005-10-30 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63266328000, #    utc_start 2005-10-30 04:00:00 (Sun)
63279637200, #      utc_end 2006-04-02 05:00:00 (Sun)
63266317200, #  local_start 2005-10-30 01:00:00 (Sun)
63279626400, #    local_end 2006-04-02 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63279637200, #    utc_start 2006-04-02 05:00:00 (Sun)
63297777600, #      utc_end 2006-10-29 04:00:00 (Sun)
63279630000, #  local_start 2006-04-02 03:00:00 (Sun)
63297770400, #    local_end 2006-10-29 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63297777600, #    utc_start 2006-10-29 04:00:00 (Sun)
63309272400, #      utc_end 2007-03-11 05:00:00 (Sun)
63297766800, #  local_start 2006-10-29 01:00:00 (Sun)
63309261600, #    local_end 2007-03-11 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63309272400, #    utc_start 2007-03-11 05:00:00 (Sun)
63329832000, #      utc_end 2007-11-04 04:00:00 (Sun)
63309265200, #  local_start 2007-03-11 03:00:00 (Sun)
63329824800, #    local_end 2007-11-04 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63329832000, #    utc_start 2007-11-04 04:00:00 (Sun)
63340722000, #      utc_end 2008-03-09 05:00:00 (Sun)
63329821200, #  local_start 2007-11-04 01:00:00 (Sun)
63340711200, #    local_end 2008-03-09 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63340722000, #    utc_start 2008-03-09 05:00:00 (Sun)
63361281600, #      utc_end 2008-11-02 04:00:00 (Sun)
63340714800, #  local_start 2008-03-09 03:00:00 (Sun)
63361274400, #    local_end 2008-11-02 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63361281600, #    utc_start 2008-11-02 04:00:00 (Sun)
63372171600, #      utc_end 2009-03-08 05:00:00 (Sun)
63361270800, #  local_start 2008-11-02 01:00:00 (Sun)
63372160800, #    local_end 2009-03-08 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63372171600, #    utc_start 2009-03-08 05:00:00 (Sun)
63392731200, #      utc_end 2009-11-01 04:00:00 (Sun)
63372164400, #  local_start 2009-03-08 03:00:00 (Sun)
63392724000, #    local_end 2009-11-01 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63392731200, #    utc_start 2009-11-01 04:00:00 (Sun)
63404226000, #      utc_end 2010-03-14 05:00:00 (Sun)
63392720400, #  local_start 2009-11-01 01:00:00 (Sun)
63404215200, #    local_end 2010-03-14 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63404226000, #    utc_start 2010-03-14 05:00:00 (Sun)
63424785600, #      utc_end 2010-11-07 04:00:00 (Sun)
63404218800, #  local_start 2010-03-14 03:00:00 (Sun)
63424778400, #    local_end 2010-11-07 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63424785600, #    utc_start 2010-11-07 04:00:00 (Sun)
63435675600, #      utc_end 2011-03-13 05:00:00 (Sun)
63424774800, #  local_start 2010-11-07 01:00:00 (Sun)
63435664800, #    local_end 2011-03-13 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63435675600, #    utc_start 2011-03-13 05:00:00 (Sun)
63456235200, #      utc_end 2011-11-06 04:00:00 (Sun)
63435668400, #  local_start 2011-03-13 03:00:00 (Sun)
63456228000, #    local_end 2011-11-06 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63456235200, #    utc_start 2011-11-06 04:00:00 (Sun)
63467125200, #      utc_end 2012-03-11 05:00:00 (Sun)
63456224400, #  local_start 2011-11-06 01:00:00 (Sun)
63467114400, #    local_end 2012-03-11 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63467125200, #    utc_start 2012-03-11 05:00:00 (Sun)
63487684800, #      utc_end 2012-11-04 04:00:00 (Sun)
63467118000, #  local_start 2012-03-11 03:00:00 (Sun)
63487677600, #    local_end 2012-11-04 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63487684800, #    utc_start 2012-11-04 04:00:00 (Sun)
63498574800, #      utc_end 2013-03-10 05:00:00 (Sun)
63487674000, #  local_start 2012-11-04 01:00:00 (Sun)
63498564000, #    local_end 2013-03-10 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63498574800, #    utc_start 2013-03-10 05:00:00 (Sun)
63519134400, #      utc_end 2013-11-03 04:00:00 (Sun)
63498567600, #  local_start 2013-03-10 03:00:00 (Sun)
63519127200, #    local_end 2013-11-03 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63519134400, #    utc_start 2013-11-03 04:00:00 (Sun)
63530024400, #      utc_end 2014-03-09 05:00:00 (Sun)
63519123600, #  local_start 2013-11-03 01:00:00 (Sun)
63530013600, #    local_end 2014-03-09 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63530024400, #    utc_start 2014-03-09 05:00:00 (Sun)
63550584000, #      utc_end 2014-11-02 04:00:00 (Sun)
63530017200, #  local_start 2014-03-09 03:00:00 (Sun)
63550576800, #    local_end 2014-11-02 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63550584000, #    utc_start 2014-11-02 04:00:00 (Sun)
63561474000, #      utc_end 2015-03-08 05:00:00 (Sun)
63550573200, #  local_start 2014-11-02 01:00:00 (Sun)
63561463200, #    local_end 2015-03-08 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63561474000, #    utc_start 2015-03-08 05:00:00 (Sun)
63582033600, #      utc_end 2015-11-01 04:00:00 (Sun)
63561466800, #  local_start 2015-03-08 03:00:00 (Sun)
63582026400, #    local_end 2015-11-01 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63582033600, #    utc_start 2015-11-01 04:00:00 (Sun)
63593528400, #      utc_end 2016-03-13 05:00:00 (Sun)
63582022800, #  local_start 2015-11-01 01:00:00 (Sun)
63593517600, #    local_end 2016-03-13 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63593528400, #    utc_start 2016-03-13 05:00:00 (Sun)
63614088000, #      utc_end 2016-11-06 04:00:00 (Sun)
63593521200, #  local_start 2016-03-13 03:00:00 (Sun)
63614080800, #    local_end 2016-11-06 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63614088000, #    utc_start 2016-11-06 04:00:00 (Sun)
63624978000, #      utc_end 2017-03-12 05:00:00 (Sun)
63614077200, #  local_start 2016-11-06 01:00:00 (Sun)
63624967200, #    local_end 2017-03-12 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63624978000, #    utc_start 2017-03-12 05:00:00 (Sun)
63645537600, #      utc_end 2017-11-05 04:00:00 (Sun)
63624970800, #  local_start 2017-03-12 03:00:00 (Sun)
63645530400, #    local_end 2017-11-05 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63645537600, #    utc_start 2017-11-05 04:00:00 (Sun)
63656427600, #      utc_end 2018-03-11 05:00:00 (Sun)
63645526800, #  local_start 2017-11-05 01:00:00 (Sun)
63656416800, #    local_end 2018-03-11 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63656427600, #    utc_start 2018-03-11 05:00:00 (Sun)
63676987200, #      utc_end 2018-11-04 04:00:00 (Sun)
63656420400, #  local_start 2018-03-11 03:00:00 (Sun)
63676980000, #    local_end 2018-11-04 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63676987200, #    utc_start 2018-11-04 04:00:00 (Sun)
63687877200, #      utc_end 2019-03-10 05:00:00 (Sun)
63676976400, #  local_start 2018-11-04 01:00:00 (Sun)
63687866400, #    local_end 2019-03-10 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63687877200, #    utc_start 2019-03-10 05:00:00 (Sun)
63708436800, #      utc_end 2019-11-03 04:00:00 (Sun)
63687870000, #  local_start 2019-03-10 03:00:00 (Sun)
63708429600, #    local_end 2019-11-03 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63708436800, #    utc_start 2019-11-03 04:00:00 (Sun)
63719326800, #      utc_end 2020-03-08 05:00:00 (Sun)
63708426000, #  local_start 2019-11-03 01:00:00 (Sun)
63719316000, #    local_end 2020-03-08 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63719326800, #    utc_start 2020-03-08 05:00:00 (Sun)
63739886400, #      utc_end 2020-11-01 04:00:00 (Sun)
63719319600, #  local_start 2020-03-08 03:00:00 (Sun)
63739879200, #    local_end 2020-11-01 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63739886400, #    utc_start 2020-11-01 04:00:00 (Sun)
63751381200, #      utc_end 2021-03-14 05:00:00 (Sun)
63739875600, #  local_start 2020-11-01 01:00:00 (Sun)
63751370400, #    local_end 2021-03-14 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63751381200, #    utc_start 2021-03-14 05:00:00 (Sun)
63771940800, #      utc_end 2021-11-07 04:00:00 (Sun)
63751374000, #  local_start 2021-03-14 03:00:00 (Sun)
63771933600, #    local_end 2021-11-07 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63771940800, #    utc_start 2021-11-07 04:00:00 (Sun)
63782830800, #      utc_end 2022-03-13 05:00:00 (Sun)
63771930000, #  local_start 2021-11-07 01:00:00 (Sun)
63782820000, #    local_end 2022-03-13 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63782830800, #    utc_start 2022-03-13 05:00:00 (Sun)
63803390400, #      utc_end 2022-11-06 04:00:00 (Sun)
63782823600, #  local_start 2022-03-13 03:00:00 (Sun)
63803383200, #    local_end 2022-11-06 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63803390400, #    utc_start 2022-11-06 04:00:00 (Sun)
63814280400, #      utc_end 2023-03-12 05:00:00 (Sun)
63803379600, #  local_start 2022-11-06 01:00:00 (Sun)
63814269600, #    local_end 2023-03-12 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63814280400, #    utc_start 2023-03-12 05:00:00 (Sun)
63834840000, #      utc_end 2023-11-05 04:00:00 (Sun)
63814273200, #  local_start 2023-03-12 03:00:00 (Sun)
63834832800, #    local_end 2023-11-05 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
    [
63834840000, #    utc_start 2023-11-05 04:00:00 (Sun)
63845730000, #      utc_end 2024-03-10 05:00:00 (Sun)
63834829200, #  local_start 2023-11-05 01:00:00 (Sun)
63845719200, #    local_end 2024-03-10 02:00:00 (Sun)
-10800,
0,
'PMST',
    ],
    [
63845730000, #    utc_start 2024-03-10 05:00:00 (Sun)
63866289600, #      utc_end 2024-11-03 04:00:00 (Sun)
63845722800, #  local_start 2024-03-10 03:00:00 (Sun)
63866282400, #    local_end 2024-11-03 02:00:00 (Sun)
-7200,
1,
'PMDT',
    ],
];

sub olson_version { '2013d' }

sub has_dst_changes { 38 }

sub _max_year { 2023 }

sub _new_instance
{
    return shift->_init( @_, spans => $spans );
}

sub _last_offset { -10800 }

my $last_observance = bless( {
  'format' => 'PM%sT',
  'gmtoff' => '-3:00',
  'local_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 725372,
    'local_rd_secs' => 0,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 725372,
    'utc_rd_secs' => 0,
    'utc_year' => 1988
  }, 'DateTime' ),
  'offset_from_std' => 0,
  'offset_from_utc' => -10800,
  'until' => [],
  'utc_start_datetime' => bless( {
    'formatter' => undef,
    'local_rd_days' => 725372,
    'local_rd_secs' => 10800,
    'offset_modifier' => 0,
    'rd_nanosecs' => 0,
    'tz' => bless( {
      'name' => 'floating',
      'offset' => 0
    }, 'DateTime::TimeZone::Floating' ),
    'utc_rd_days' => 725372,
    'utc_rd_secs' => 10800,
    'utc_year' => 1988
  }, 'DateTime' )
}, 'DateTime::TimeZone::OlsonDB::Observance' )
;
sub _last_observance { $last_observance }

my $rules = [
  bless( {
    'at' => '2:00',
    'from' => '2007',
    'in' => 'Nov',
    'letter' => 'S',
    'name' => 'Canada',
    'offset_from_std' => 0,
    'on' => 'Sun>=1',
    'save' => '0',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' ),
  bless( {
    'at' => '2:00',
    'from' => '2007',
    'in' => 'Mar',
    'letter' => 'D',
    'name' => 'Canada',
    'offset_from_std' => 3600,
    'on' => 'Sun>=8',
    'save' => '1:00',
    'to' => 'max',
    'type' => undef
  }, 'DateTime::TimeZone::OlsonDB::Rule' )
]
;
sub _rules { $rules }


1;

