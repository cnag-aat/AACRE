package PDF::API2::Page;

use base 'PDF::API2::Basic::PDF::Pages';

use strict;
use warnings;

our $VERSION = '2.035'; # VERSION

use POSIX qw(floor);
use Scalar::Util qw(weaken);

use PDF::API2::Basic::PDF::Utils;
use PDF::API2::Content;
use PDF::API2::Content::Text;
use PDF::API2::Util;

=head1 NAME

PDF::API2::Page - Methods to interact with individual pages

=head1 METHODS

=over

=item $page = PDF::API2::Page->new $pdf, $parent, $index

Returns a page object (called from $pdf->page).

=cut

sub new {
    my ($class, $pdf, $parent, $index) = @_;
    my $self = {};

    $class = ref($class) if ref($class);
    $self = $class->SUPER::new($pdf, $parent);
    $self->{'Type'} = PDFName('Page');
    $self->proc_set(qw( PDF Text ImageB ImageC ImageI ));
    delete $self->{'Count'};
    delete $self->{'Kids'};
    $parent->add_page($self, $index);
    return $self;
}

=item $page = PDF::API2::Page->coerce $pdf, $pdfpage

Returns a page object converted from $pdfpage (called from $pdf->openpage).

=cut

sub coerce {
    my ($class, $pdf, $page) = @_;
    my $self = $page;
    bless $self, $class;
    $self->{' apipdf'} = $pdf;
    weaken $self->{' apipdf'};
    return $self;
}

=item $page->update

Marks a page to be updated (by $pdf->update).

=cut

sub update {
    my ($self) = @_;
    $self->{' apipdf'}->out_obj($self);
    return $self;
}

=item ($llx, $lly, $urx, $ury) = $page->mediabox()

=item $page->mediabox($w, $h)

=item $page->mediabox($llx, $lly, $urx, $ury)

=item $page->mediabox($alias)

Get or set the mediabox.  This method supports the following aliases:
'4A0', '2A0', 'A0', 'A1', 'A2', 'A3', 'A4', 'A5', 'A6',
'4B0', '2B0', 'B0', 'B1', 'B2', 'B3', 'B4', 'B5', 'B6',
'LETTER', 'BROADSHEET', 'LEDGER', 'TABLOID', 'LEGAL',
'EXECUTIVE', and '36X36'.

=cut

sub _bounding_box {
    my $self = shift();
    my $type = shift();

    # Get
    unless (scalar @_) {
        my $box = $self->find_prop($type);
        unless ($box) {
            # Default to letter (for historical PDF::API2 reasons, not per the
            # PDF specification)
            return (0, 0, 612, 792) if $type eq 'MediaBox';

            # Use defaults per PDF 1.7 section 14.11.2 Page Boundaries
            return $self->_bounding_box('MediaBox') if $type eq 'CropBox';
            return $self->_bounding_box('CropBox');
        }
        return map { $_->val() } $box->elements();
    }

    # Set
    $self->{$type} = PDFArray(map { PDFNum(float($_)) } page_size(@_));
    return $self;
}

sub mediabox {
    my $self = shift();
    return $self->_bounding_box('MediaBox', @_);
}

# Deprecated
sub get_mediabox {
    my $self = shift();
    return $self->_bounding_box('MediaBox');
}

=item ($llx, $lly, $urx, $ury) = $page->cropbox()

=item $page->cropbox($w, $h)

=item $page->cropbox($llx, $lly, $urx, $ury)

=item $page->cropbox($alias)

Get or set the cropbox.  This method supports the same aliases as mediabox.

The cropbox defaults to the mediabox.

=cut

sub cropbox {
    my $self = shift();
    return $self->_bounding_box('CropBox', @_);
}

# Deprecated
sub get_cropbox {
    my $self = shift();
    return $self->_bounding_box('CropBox');
}

=item ($llx, $lly, $urx, $ury) = $page->bleedbox()

=item $page->bleedbox($w, $h)

=item $page->bleedbox($llx, $lly, $urx, $ury)

=item $page->bleedbox($alias)

Get or set the bleedbox.  This method supports the same aliases as mediabox.

The bleedbox defaults to the cropbox.

=cut

sub bleedbox {
    my $self = shift();
    return $self->_bounding_box('BleedBox', @_);
}

# Deprecated
sub get_bleedbox {
    my $self = shift();
    return $self->_bounding_box('BleedBox');
}

=item ($llx, $lly, $urx, $ury) = $page->trimbox()

=item $page->trimbox($w, $h)

=item $page->trimbox($llx, $lly, $urx, $ury)

Get or set the trimbox.  This method supports the same aliases as mediabox.

The trimbox defaults to the cropbox.

=cut

sub trimbox {
    my $self = shift();
    return $self->_bounding_box('TrimBox', @_);
}

# Deprecated
sub get_trimbox {
    my $self = shift();
    return $self->_bounding_box('TrimBox');
}

=item ($llx, $lly, $urx, $ury) = $page->artbox()

=item $page->artbox($w, $h)

=item $page->artbox($llx, $lly, $urx, $ury)

=item $page->artbox($alias)

Get or set the artbox.  This method supports the same aliases as mediabox.

The rtbox defaults to the cropbox.

=cut

sub artbox {
    my $self = shift();
    return $self->_bounding_box('ArtBox', @_);
}

# Deprecated
sub get_artbox {
    my $self = shift();
    return $self->_bounding_box('ArtBox');
}

=item $page->rotate $deg

Rotates the page by the given degrees, which must be a multiple of 90.

(This allows you to auto-rotate to landscape without changing the mediabox!)

=cut

sub rotate {
    my ($self, $degrees) = @_;

    # Ignore rotation of 360 or more (in either direction)
    $degrees = $degrees % 360;

    $self->{'Rotate'} = PDFNum($degrees);

    return $self;
}

=item $gfx = $page->gfx $prepend

Returns a graphics content object. If $prepend is true the content
will be prepended to the page description.

=cut

sub fixcontents {
    my $self = shift();
    $self->{'Contents'} = $self->{'Contents'} || PDFArray();
    if (ref($self->{'Contents'}) =~ /Objind$/) {
        $self->{'Contents'}->realise();
    }
    if (ref($self->{'Contents'}) !~ /Array$/) {
        $self->{'Contents'} = PDFArray($self->{'Contents'});
    }
    return;
}

sub content {
    my ($self, $obj, $dir) = @_;
    if (defined($dir) and $dir > 0) {
        $self->precontent($obj);
    }
    else {
        $self->addcontent($obj);
    }
    $self->{' apipdf'}->new_obj($obj) unless $obj->is_obj($self->{' apipdf'});
    $obj->{' apipdf'} = $self->{' apipdf'};
    $obj->{' api'} = $self->{' api'};
    $obj->{' apipage'} = $self;

    weaken $obj->{' apipdf'};
    weaken $obj->{' api'};
    weaken $obj->{' apipage'};

    return $obj;
}

sub addcontent {
    my ($self, @objs) = @_;
    $self->fixcontents();
    $self->{'Contents'}->add_elements(@objs);
    return;
}
sub precontent {
    my ($self, @objs) = @_;
    $self->fixcontents();
    unshift @{$self->{'Contents'}->val()}, @objs;
    return;
}

sub gfx {
    my ($self, $dir) = @_;
    my $gfx=PDF::API2::Content->new();
    $self->content($gfx,$dir);
    $gfx->compressFlate() if $self->{' api'}->{'forcecompress'};
    return $gfx;
}

=item $txt = $page->text $prepend

Returns a text content object. If $prepend is true the content
will be prepended to the page description.

=cut

sub text {
    my ($self, $dir) = @_;
    my $text = PDF::API2::Content::Text->new();
    $self->content($text, $dir);
    $text->compressFlate() if $self->{' api'}->{'forcecompress'};
    return $text;
}

=item $ant = $page->annotation

Returns a new annotation object.

=cut

sub annotation {
    my $self = shift();

    unless (exists $self->{'Annots'}) {
        $self->{'Annots'} = PDFArray();
        $self->update();
    }
    elsif (ref($self->{'Annots'}) =~ /Objind/) {
        $self->{'Annots'}->realise();
    }

    require PDF::API2::Annotation;
    my $ant = PDF::API2::Annotation->new();
    $self->{'Annots'}->add_elements($ant);
    $self->{' apipdf'}->new_obj($ant);
    $ant->{' apipdf'} = $self->{' apipdf'};
    $ant->{' apipage'} = $self;
    weaken $ant->{' apipdf'};
    weaken $ant->{' apipage'};

    if ($self->{'Annots'}->is_obj($self->{' apipdf'})) {
        $self->{' apipdf'}->out_obj($self->{'Annots'});
    }

    return $ant;
}

=item $page->resource $type, $key, $obj

Adds a resource to the page-inheritance tree.

B<Example:>

    $co->resource('Font',$fontkey,$fontobj);
    $co->resource('XObject',$imagekey,$imageobj);
    $co->resource('Shading',$shadekey,$shadeobj);
    $co->resource('ColorSpace',$spacekey,$speceobj);

B<Note:> You only have to add the required resources, if
they are NOT handled by the *font*, *image*, *shade* or *space*
methods.

=cut

sub resource {
    my ($self, $type, $key, $obj, $force) = @_;
    my $dict = $self->find_prop('Resources');

    $dict = $dict || $self->{'Resources'} || PDFDict();

    $dict->realise() if ref($dict) =~ /Objind$/;

    $dict->{$type} = $dict->{$type} || PDFDict();
    $dict->{$type}->realise if ref($dict->{$type}) =~ /Objind$/;

    unless (defined $obj) {
        return $dict->{$type}->{$key} || undef;
    }
    else {
        if ($force) {
            $dict->{$type}->{$key} = $obj;
        }
        else {
            $dict->{$type}->{$key} = $dict->{$type}->{$key} || $obj;
        }

        $self->{' apipdf'}->out_obj($dict)          if $dict->is_obj($self->{' apipdf'});
        $self->{' apipdf'}->out_obj($dict->{$type}) if $dict->{$type}->is_obj($self->{' apipdf'});
        $self->{' apipdf'}->out_obj($obj)           if $obj->is_obj($self->{' apipdf'});
        $self->{' apipdf'}->out_obj($self);

        return $dict;
    }
}

sub ship_out {
    my ($self, $pdf) = @_;

    $pdf->ship_out($self);
    if (defined $self->{'Contents'}) {
        $pdf->ship_out($self->{'Contents'}->elements());
    }
    return $self;
}

=back

=cut

1;
