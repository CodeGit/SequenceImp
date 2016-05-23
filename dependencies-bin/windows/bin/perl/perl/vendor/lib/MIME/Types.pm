# Copyrights 1999,2001-2013 by [Mark Overmeer].
#  For other contributors see ChangeLog.
# See the manual pages for details on the licensing terms.
# Pod stripped from pm file by OODoc 2.01.

package MIME::Types;
use vars '$VERSION';
$VERSION = '2.01';


use strict;

use MIME::Type     ();
use File::Spec     ();
use File::Basename qw(dirname);


my %typedb;
sub new(@) { (bless {}, shift)->init( {@_} ) }

sub init($)
{   my ($self, $args) = @_;
    keys %typedb or $self->_read_db($args);
    $self;
}

sub _read_db($)
{   my ($self, $args)   = @_;
    my $skip_extensions = $args->{skip_extensions};
    my $only_complete   = $args->{only_complete};
    my $only_iana       = $args->{only_iana};

    my $db              = $args->{db_file}
      || File::Spec->catfile(dirname(__FILE__), 'types.db');

    open DB, '<:encoding(utf8)', $db
       or die "cannot open type database in $db: $!\n";

    while(1)
    {   my $header = <DB>;
        defined $header or last;
        chomp $header;

        # This logic is entangled with the bin/collect_types script
        my ($count, $major, $is_iana, $has_ext) = split /\:/, $header;
        my $skip_section = $major eq 'EXTENSIONS' ? $skip_extensions
          : (($only_iana && !$is_iana) || ($only_complete && !$has_ext));

#warn "Skipping section $header\n" if $skip_section;
        (my $section = $major) =~ s/^x-//;
        if($major eq 'EXTENSIONS')
        {   while(<DB>)
            {   last if m/^$/;
                next if $skip_section;
                chomp;
                $typedb{$section}{$1} = $2 if m/(.*);(.*)/;
            }
        }
        else
        {   while(<DB>)
            {   last if m/^$/;
                next if $skip_section;
                chomp;
                $typedb{$section}{$1} = "$major/$_" if m/^(?:x-)?([^;]+)/;
            }
        }
    }

    close DB;
}

# Catalyst-Plugin-Static-Simple uses it :(
sub create_type_index {}

#-------------------------------------------


sub type($)
{   my $spec    = lc $_[1];
    $spec       = 'text/plain' if $spec eq 'text';   # old mailers

    $spec =~ m!^(?:x\-)?([^/]+)/(?:x-)?(.*)!
        or return;

    my $section = $typedb{$1}    or return;
    my $record  = $section->{$2} or return;
    return $record if ref $record;   # already extended

    my $simple   = $2;
    my ($type, $ext, $enc) = split m/\;/, $record;
    my $os       = undef;   # XXX TODO

    $section->{$simple} = MIME::Type->new
      ( type       => $type
      , extensions => [split /\,/, $ext]
      , encoding   => $enc
      , system     => $os
      );
}


sub mimeTypeOf($)
{   my ($self, $name) = @_;
    (my $ext = lc $name) =~ s/.*\.//;
    my $type = $typedb{EXTENSIONS}{$ext} or return;
    $self->type($type);
}


sub addType(@)
{   my $self = shift;

    foreach my $type (@_)
    {   my ($major, $minor) = split m!/!, $type->simplified;
        $typedb{$major}{$minor} = $type;
        $typedb{EXTENSIONS}{$_} = $type for $type->extensions;
    }
    $self;
}


sub types()
{   my $self  = shift;
    my @types;
    foreach my $section (keys %typedb)
    {   next if $section eq 'EXTENSIONS';
        push @types, map $_->type("$section/$_"),
                         sort keys %{$typedb{$section}};
    }
    @types;
}


sub listTypes()
{   my $self  = shift;
    my @types;
    foreach my $section (keys %typedb)
    {   next if $section eq 'EXTENSIONS';
        foreach my $sub (sort keys %{$typedb{$section}})
        {   my $record = $typedb{$section}{$sub};
            push @types, ref $record            ? $record->type
                       : $record =~ m/^([^;]+)/ ? $1 : die;
        }
    }
    @types;
}



sub extensions { keys %{$typedb{EXTENSIONS}} }

#-------------------------------------------
# OLD INTERGFACE (version 0.06 and lower)


use base 'Exporter';
our @EXPORT_OK = qw(by_suffix by_mediatype import_mime_types);


my $mime_types;

sub by_suffix($)
{   my $filename = shift;
    $mime_types ||= MIME::Types->new;
    my $mime     = $mime_types->mimeTypeOf($filename);

    my @data     = defined $mime ? ($mime->type, $mime->encoding) : ('','');
    wantarray ? @data : \@data;
}


sub by_mediatype($)
{   my $type = shift;
    $mime_types ||= MIME::Types->new;

    my @found;
    if(!ref $type && index($type, '/') >= 0)
    {   my $mime   = $mime_types->type($type);
        @found     = $mime if $mime;
    }
    else
    {   my $search = ref $type eq 'Regexp' ? $type : qr/$type/i;
        @found     = map $mime_types->type($_),
                         grep $_ =~ $search,
                             $mime_types->listTypes;
    }

    my @data;
    foreach my $mime (@found)
    {   push @data, map [$_, $mime->type, $mime->encoding],
                        $mime->extensions;
    }

    wantarray ? @data : \@data;
}


sub import_mime_types($)
{   my $filename = shift;
    use Carp;
    croak <<'CROAK';
import_mime_types is not supported anymore: if you have types to add
please send them to the author.
CROAK
}

1;
__END__
# Exceptions
vms:text/plain;doc;8bit
mac:application/x-macbase64;;bin

# IE6 bug
image/pjpeg;;base64
