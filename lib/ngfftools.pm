package ngfftools;
#no warnings 'experimental::smartmatch';
use Data::Dumper;
use strict;
use Time::HiRes qw/time/;
use FindBin qw($RealBin);
use Exporter;
use Getopt::Long;
use Cwd qw/abs_path/;
our @ISA       = qw(Exporter);
our @EXPORT_OK = qw($parsers $autofillers $printers detector_informat detector_outformat isPossibleCombination merger abs2rel rel2abs seq extract);

my $json;
my $Levels;
my %Gene;

sub ctime {
    chomp( my $tmp = `date` );
    return $tmp;
}

sub timer {
    my $prefix = shift;
    $prefix ||= "";
    chomp( my $tmp = `date` );
    print STDERR $prefix . " : " . $tmp . "\n";
}

sub extract {
    my ( $json, $UNIT, $outfh ) = @_;
    my %nameChilds = (
        'chromosome' => 'locus',
        'locus'      => 'primary',
        'primary'    => 'processed',
        'processed'  => 'product'
    );
    my %nameParents = reverse %nameChilds;
    my %Done;
    my %Record;
    my @Level = ( sort keys %nameChilds, "product" );
    foreach my $level (@Level) {
        foreach my $id ( sort keys %{$UNIT} ) {
            if ( exists $json->{$level}->{ids}->{$id} ) {
                my @str = ( $id, $level, $outfh );
                push @{ $Record{$id} }, \@str;
                if ( $UNIT->{$id} eq "parent" ) {
                    my @parentIds   = sort keys %{ $json->{$level}->{ids}->{$id}->{parentId} };
                    my $parentLevel = $nameParents{$level};
                    for my $parentId (@parentIds) {
                        my @str = ( $parentId, $parentLevel, $outfh );
                        push @{ $Record{$id} }, \@str;
                        if ( exists $nameParents{$parentLevel} ) {
                            my @parentIds2   = sort keys %{ $json->{$parentLevel}->{ids}->{$parentId}->{parentId} };
                            my $parentLevel2 = $nameParents{$parentLevel};
                            for my $parentId2 (@parentIds2) {
                                my @str = ( $parentId2, $parentLevel2, $outfh );
                                push @{ $Record{$id} }, \@str;
                                if ( exists $nameParents{$parentLevel2} ) {
                                    my @parentIds3   = sort keys %{ $json->{$parentLevel2}->{ids}->{$parentId2}->{parentId} };
                                    my $parentLevel3 = $nameParents{$parentLevel2};
                                    for my $parentId3 (@parentIds3) {
                                        my @str = ( $parentId3, $parentLevel3, $outfh );
                                        push @{ $Record{$id} }, \@str;
                                        if ( exists $nameParents{$parentLevel3} ) {
                                            my @parentIds4   = sort keys %{ $json->{$parentLevel3}->{ids}->{$parentId3}->{parentId} };
                                            my $parentLevel4 = $nameParents{$parentLevel3};
                                            for my $parentId4 (@parentIds4) {
                                                my @str = ( $parentId4, $parentLevel4, $outfh );
                                                push @{ $Record{$id} }, \@str;
                                            }
                                        }
                                    }
                                }
                            }
                        }
                    }
                }
                if ( $UNIT->{$id} eq "child" ) {
                    my @childIds   = sort keys %{ $json->{$level}->{ids}->{$id}->{childId} };
                    my $childLevel = $nameChilds{$level};
                    for my $childId (@childIds) {
                        my @str = ( $childId, $childLevel, $outfh );
                        push @{ $Record{$id} }, \@str;
                        if ( exists $nameChilds{$childLevel} ) {
                            my @childIds2   = sort keys %{ $json->{$childLevel}->{ids}->{$childId}->{childId} };
                            my $childLevel2 = $nameChilds{$childLevel};
                            for my $childId2 (@childIds2) {
                                my @str = ( $childId2, $childLevel2, $outfh );
                                push @{ $Record{$id} }, \@str;
                                if ( exists $nameChilds{$childLevel2} ) {
                                    my @childIds3   = sort keys %{ $json->{$childLevel2}->{ids}->{$childId2}->{childId} };
                                    my $childLevel3 = $nameChilds{$childLevel2};
                                    for my $childId3 (@childIds3) {
                                        my @str = ( $childId3, $childLevel3, $outfh );
                                        push @{ $Record{$id} }, \@str;
                                        if ( exists $nameChilds{$childLevel3} ) {
                                            my @childIds4   = sort keys %{ $json->{$childLevel3}->{ids}->{$childId3}->{childId} };
                                            my $childLevel4 = $nameChilds{$childLevel3};
                                            for my $childId4 (@childIds4) {
                                                my @str = ( $childId4, $childLevel4, $outfh );
                                                push @{ $Record{$id} }, \@str;
                                            }
                                        }
                                    }
                                }

                            }
                        }
                    }

                }
            }
        }
    }

    foreach my $id ( sort keys %Record ) {
        my @values = @{ $Record{$id} };
        if ( $UNIT->{$id} eq "parent" ) {
            @values = reverse(@values);
        }
        for my $each (@values) {
            my @str = @$each;
            printer_subngff(@str);
        }
    }
}

sub rel2abs {
    my ( $RELposFile, $outfh ) = @_;
    open REL, "<$RELposFile";
    while (<REL>) {
        chomp;
        next if /^#/;
        my ( $seqId, $start, $end, $strand, $level, $featureType, $id, $parentId, $attributes ) = split /\t/;
        my @outpos = retrieve( $id, $level, abs($start), abs($end), "rel2abs", $json );
        for my $pos (@outpos) {
            my ( $absstart, $absend ) = split /\s/, $pos;
            next if !defined $absstart || !defined $absend;
            if ( $absstart > $absend ) {
                my $tmp = $absstart;
                $absstart = $absend;
                $absend   = $tmp;
            }
            print $outfh join( "\t", $seqId, $absstart, $absend, $strand, $level, $featureType, $id, $parentId, $attributes ) . "\n";
        }
    }
    close REL;
}

sub seq {
    my ( $json, $genomeFile, $species_type, $outpre ) = @_;
    my %Genome;
    if ($genomeFile =~ /\.gz$/i){
        open GENOME, "gzip -dc $genomeFile|" or die "Cannot open file $genomeFile: $!\n";
    }else{
        open GENOME, "<$genomeFile" or die "Cannot open file $genomeFile: $!\n";
    }
    $/ = ">";
    <GENOME>;
    while (<GENOME>) {
        chomp;
        my ( $header, $seqs ) = split /\n/, $_, 2;
        $seqs   =~ s/\s+//sg;
        $header =~ s/^(\S+) .*/$1/;
        $Genome{$header} = $seqs;
    }
    $/ = "\n";
    close GENOME;
    my @levels_order = ( "chromosome", "locus", "primary", "processed", "product" );
    my ( $OUTCDS, $OUTPEP );
    open OUTNUC, ">${outpre}_nucleotide.fasta";
    open $OUTCDS, ">${outpre}_CDS.fasta";
    open $OUTPEP, ">${outpre}_pep.fasta";
    open OUTREG, ">${outpre}_regulator.fasta";
    my %Chr_type;
    my %nameTrans = (
        'major'                      => 'id1',
        'Chloroplast'                => 'id11',
        'Vertebrate_mitochondrion'   => 'id2',
        'Yeast_mitochondrion'        => 'id3',
        'Protozoan_mitochondrion'    => 'id4',
        'Invertebrate_mitochondrion' => 'id5'
    );

    for my $level (@levels_order) {
        for my $id ( @{ $json->{$level}->{turns} } ) {
            if ( exists $json->{$level}->{"ids"}->{$id}->{nonRegulator} ) {
                my @Values    = @{ $json->{$level}->{"ids"}->{$id}->{nonRegulator} };
                my @ValuesCor = @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulatorCor"} };
                if ( $level eq "chromosome" ) {
                    my $attr = @{ $Values[0] }[5];
                    $attr =~ s/seqSource=([^;]*?);.*/$1/;
                    $Chr_type{$id} = $attr;
                }
                else {
                    my ( $seqs, $chr, $start, $end, $strand, $type,$editing);
		    if (exists $json->{$level}->{ids}->{$id}->{nonRegulator_shared_attributes}->{edited_site}){
			$editing=$json->{$level}->{ids}->{$id}->{nonRegulator_shared_attributes}->{edited_site};
		    }
                    my %hash_tmp;
                    foreach my $i ( 0 .. $#ValuesCor ) {
                        my $each = $ValuesCor[$i];
                        my ( $rstart, $rend ) = @$each;
                        $hash_tmp{$rstart}{indx} = $i;
                    }
                    foreach my $rstart ( sort { $a <=> $b } keys %hash_tmp ) {
                        my $indx = $hash_tmp{$rstart}{indx};

                        ( $chr, $start, $end, $strand, $type ) = @{ $Values[$indx] }[ 0, 1, 2, 3, 4];
                        my $seq = uc( substr( $Genome{$chr}, $start - 1, $end - $start + 1 ) );
                        if ( $strand eq "-" ) {

                            #because minus strand has been ordered ascendingly in nonRegulatorCor, so do not need reverse.
                            $seq =~ tr/ATGC/TACG/;
                        }
			if(defined $editing){
			    my @edited_sites=split /,/,$editing;
			    foreach my $site (@edited_sites){
				my ($editingChr,$editingstring)=split /:/,$site;
				if ($editingstring=~/[A-Z](\d+)([A-Z])/){
					my ($editingPos,$editingNucl)=($1,$2);
					next if $editingPos < $start || $editingPos > $end;
			                $seq=edited_site($start,$end,$strand,$seq,$editingPos,$editingNucl);
				}
			    }
			}
                        $seqs .= $seq;
                    }

		    sub edited_site{
			my ($start,$end,$strand,$seq,$editingPos,$editingNucl)=@_;
			my $relPos;
			if ($strand eq "+"){
				$relPos=$editingPos-$start;
				substr($seq,$relPos,1,$editingNucl);
			}else{
				$relPos=$end-$editingPos;
				$editingNucl=~ tr/ATGC/TACG/;
				substr($seq,$relPos,1,$editingNucl);
			}
				
			return $seq;
		    }

                    next if $start < 0;
                    if ( $type eq "CDS" ) {
                        my $code_type;
                        if ( $Chr_type{$chr} eq "major" ) {
                            $code_type = $nameTrans{"major"};
                        }
                        elsif ( $Chr_type{$chr} eq "Chloroplast" ) {
                            $code_type = $nameTrans{"Chloroplast"};
                        }
                        elsif ( $Chr_type{$chr} eq "mitochondrion" ) {
                            $code_type = $nameTrans{"${species_type}_mitochondrion"};
                        }
                        else {
                            print STDERR "no special codon system for $chr($Chr_type{$chr}); using id1 system\n";
                            $code_type = $nameTrans{"major"};
                        }
                        if ( !exists $Genome{$chr} ) {
                            print STDERR "$chr don't have sequence in genome file!\n";
                            next;
                        }
                        print STDERR "#$chr use codon $code_type\n";
                        CDS2PEP( $seqs, $id, $code_type, $outpre, $OUTCDS, $OUTPEP );
                    }
                    if ( !exists $Genome{$chr} ) {
                        print STDERR "$chr don't have sequence in genome file!\n";
                        next;
                    }
                    print OUTNUC ">$id $level\n$seqs\n";
                }
            }
        }
    }
    close OUTNUC;
    close $OUTCDS;
    close $OUTPEP;
    for my $level (@levels_order) {
        for my $id ( @{ $json->{$level}->{turns} } ) {
            if ( exists $json->{$level}->{"ids"}->{$id}->{regulator} ) {
                foreach my $regulator_biotype ( sort keys %{ $json->{$level}->{"ids"}->{$id}->{regulator} } ) {
                    my @Pos = @{ $json->{$level}->{"ids"}->{$id}->{regulator}->{$regulator_biotype} };
                    my ( $seqs, $chr, $start, $end, $strand, $type );
                    if ( $regulator_biotype eq "polyA_sequence" ) {
                        my %hash = ArrayToHash( ${ $json->{$level}->{"ids"}->{$id}->{regulator}->{$regulator_biotype} }[0]->[5] );
                        $seqs = $hash{sequence};
                        print OUTREG ">${id}_$regulator_biotype $level\n$seqs\n";
                        next;
                    }
                    for my $i ( 0 .. $#Pos ) {
                        ( $chr, $start, $end, $strand, $type ) = @{ $Pos[$i] }[ 0, 1, 2, 3, 4 ];
                        my $seq = uc( substr( $Genome{$chr}, $start - 1, $end - $start + 1 ) );
                        if ( $strand eq "-" ) {
                            $seq = reverse($seq);
                            $seq =~ tr/ATGC/TACG/;
                        }
                        $seqs .= $seq;
                    }
                    print OUTREG ">${id}_$regulator_biotype $level\n$seqs\n";
                }
            }
        }
    }

    sub CDS2PEP {
        my %CODE = (
            'id1' => {
                'GCA' => 'A',
                'GCC' => 'A',
                'GCG' => 'A',
                'GCT' => 'A',                                                                          # Alanine
                'TGC' => 'C', 'TGT' => 'C',                                                            # Cysteine
                'GAC' => 'D', 'GAT' => 'D',                                                            # Aspartic Acid
                'GAA' => 'E', 'GAG' => 'E',                                                            # Glutamic Acid
                'TTC' => 'F', 'TTT' => 'F',                                                            # Phenylalanine
                'GGA' => 'G', 'GGC' => 'G', 'GGG' => 'G', 'GGT' => 'G',                                # Glycine
                'CAC' => 'H', 'CAT' => 'H',                                                            # Histidine
                'ATA' => 'I', 'ATC' => 'I', 'ATT' => 'I',                                              # Isoleucine
                'AAA' => 'K', 'AAG' => 'K',                                                            # Lysine
                'CTA' => 'L', 'CTC' => 'L', 'CTG' => 'L', 'CTT' => 'L', 'TTA' => 'L', 'TTG' => 'L',    # Leucine
                'ATG' => 'M',                                                                          # Methionine/Start
                'AAC' => 'N', 'AAT' => 'N',                                                            # Asparagine
                'CCA' => 'P', 'CCC' => 'P', 'CCG' => 'P', 'CCT' => 'P',                                # Proline
                'CAA' => 'Q', 'CAG' => 'Q',                                                            # Glutamine
                'CGA' => 'R', 'CGC' => 'R', 'CGG' => 'R', 'CGT' => 'R', 'AGA' => 'R', 'AGG' => 'R',    # Arginine
                'TCA' => 'S', 'TCC' => 'S', 'TCG' => 'S', 'TCT' => 'S', 'AGC' => 'S', 'AGT' => 'S',    # Serine
                'ACA' => 'T', 'ACC' => 'T', 'ACG' => 'T', 'ACT' => 'T',                                # Threonine
                'GTA' => 'V', 'GTC' => 'V', 'GTG' => 'V', 'GTT' => 'V',                                # Valine
                'TGG' => 'W',                                                                          # Tryptophan
                'TAC' => 'Y', 'TAT' => 'Y',                                                            # Tyrosine
                'TAA' => '*', 'TAG' => '*', 'TGA' => '*',                                              # Stop
                'NNN' => 'X'
            },
            'id2' => {
                'AGA' => '*',
                'AGG' => '*',
                'ATA' => 'M',
                'TGA' => 'W'
            },
            'id3' => {
                'ATA' => 'M',
                'CTT' => 'T',
                'CTC' => 'T',
                'CTA' => 'T',
                'CTG' => 'T',
                'TGA' => 'W'
            },
            'id4' => {
                'TGA' => 'W'
            },
            'id5' => {
                'AGA' => 'S',
                'AGG' => 'S',
                'ATA' => 'M',
                'TGA' => 'W'
            },
            'id11' => {}
        );

        my ( $string, $id, $code_type, $outpre, $OUTCDS, $OUTPEP ) = @_;
        my $len = length($string);
        if ( $len % 3 != 0 ) {
            print STDERR "WRONG_LENGTH\t$id\n";
        }

        sub generate {
            my %Abbrev = (
                'M' => [ 'A', 'C' ],
                'R' => [ 'A', 'G' ],
                'W' => [ 'A', 'T' ],
                'S' => [ 'C', 'G' ],
                'Y' => [ 'C', 'T' ],
                'K' => [ 'G', 'T' ],
                'V' => [ 'A', 'C', 'G' ],
                'H' => [ 'A', 'C', 'T' ],
                'D' => [ 'A', 'G', 'T' ],
                'B' => [ 'C', 'G', 'T' ],
                'X' => [ 'A', 'C', 'G', 'T' ],
                'N' => [ 'A', 'C', 'G', 'T' ]
            );
            if ( $_[0] =~ /(.*)([RYSWKBDHVN])(.*)/ ) {
                my $head = $1;
                my $tail = $3;
                my @seqs;
                foreach my $nuc ( @{ $Abbrev{$2} } ) {
                    push @seqs, generate( $head . $nuc . $tail );
                }
                return @seqs;
            }
            else {
                return $_[0];
            }
        }

        # Demo: print all sequences generated from ANCRG.
        my @nucSeq = generate($string);
        for my $i ( 0 .. $#nucSeq ) {
            my $seq = $nucSeq[$i];
            my $head;
            if ( $#nucSeq > 0 ) {
                $head = "${id}_$i";
            }
            else {
                $head = $id;
            }
            my $protein = 1;
            while (1) {
                last unless ( $seq =~ s/^(\w{3})// );
                my $sequence = uc $1;
                my $amino_a;
                if ( $sequence =~ /[^N]{3}/ ) {
                    if ( exists $CODE{$code_type}{$sequence} ) {
                        $amino_a = $CODE{$code_type}{$sequence};
                    }
                    else {
                        $amino_a = $CODE{'id1'}{$sequence};
                    }
                }
                else {
                    $amino_a = "X";
                }
                $protein .= $amino_a;
            }
            $protein =~ s/^1//;
            if ( !( $protein =~ /^M/ ) ) {
                print STDERR "warning: without start codon\t$head\n";
            }
            if ( !( $protein =~ /\*$/ ) ) {
                print STDERR "warning: without stop codon\t$head\n";
            }
            if ( $protein =~ /\*\w/ || $protein =~ /\*\*/ ) {
                print STDERR "warning: stop codon inside\t$head\n";
            }
            $protein =~ s/.{60}/$&\n/g;
            $protein =~ s/\n$//;
            print $OUTCDS ">$head\n$nucSeq[$i]\n";
            print $OUTPEP ">$head\n$protein\n";
        }
    }
}

sub retrieve {
    my ( $id, $level, $pos1, $pos2, $type, $json ) = @_;
    if ( $type eq "rel2abs" ) {
        if ( exists $json->{$level}->{"ids"}->{$id}->{nonRegulatorCor} ) {
            my ( $abspos1, $abspos2 );
            my ( $edge1, $edge2 ) = ( 0, 0 );
            my @outpos;
            my @Pos    = @{ $json->{$level}->{"ids"}->{$id}->{nonRegulatorCor} };
            my $strand = $json->{$level}->{"ids"}->{$id}->{"strand"};
	    if ( $strand eq "-" ) {
	        @Pos = reverse(@Pos);
	    }
            my @last_pos = @{ $Pos[0] };
            my ( $last_start, $last_end ) = @last_pos;
            for my $i ( 0 .. $#Pos ) {
                my ( $start, $end ) = @{ $Pos[$i] };
                my $absStart = ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[$i]->[1];
                my $absEnd   = ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[$i]->[2];
                if ( $strand eq "-" ) {
                    $absStart = ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[ $#{ $json->{$level}->{"ids"}->{$id}->{nonRegulator} } - $i ]->[1];
                    $absEnd   = ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[ $#{ $json->{$level}->{"ids"}->{$id}->{nonRegulator} } - $i ]->[2];
                }
                if (   ( ( $pos1 >= $start && $pos1 <= $end ) || ( $pos1 >= $end && $pos1 <= $start ) )
                    && ( ( $pos2 >= $start && $pos2 <= $end ) || ( $pos2 >= $end && $pos2 <= $start ) ) )
                {
                    $abspos1 = abs( $start - $pos1 ) + $absStart;
                    $abspos2 = abs( $start - $pos2 ) + $absStart;
                    push @outpos, join( " ", $abspos1, $abspos2 );
                }
                elsif ( ( $pos1 >= $start && $pos1 <= $end ) || ( $pos1 >= $end && $pos1 <= $start ) ) {
                    $abspos1 = abs( $start - $pos1 ) + $absStart;
                    if ( $end > $start ) {
                        $edge1 = $absEnd;
                        push @outpos, join( " ", $abspos1, $edge1 );
                    }
                    else {
                        $edge1 = $absStart;
                        push @outpos, join( " ", $edge1, $abspos1 );
                    }
                }
                elsif ( ( $pos2 >= $start && $pos2 <= $end ) || ( $pos2 >= $end && $pos2 <= $start ) ) {
                    $abspos2 = abs( $start - $pos2 ) + $absStart;
                    if ( $end > $start ) {
                        $edge2 = $absStart;
                        push @outpos, join( " ", $edge2, $abspos2 );

                    }
                    else {
                        $edge2 = $absEnd;
                        push @outpos, join( " ", $abspos2, $edge2 );
                    }
                }
                elsif ( $edge1 != 0 && $edge2 == 0 ) {
                    push @outpos, join( " ", $absStart, $absEnd );
                }
            }
            return @outpos;
        }
        else {
            print STDERR "$id\t$level records is wrong\n";
            return "-1";
        }
    }
    elsif ( $type eq "abs2rel" ) {
        if ( exists $json->{$level}->{"ids"}->{$id}->{nonRegulator} ) {
            my @Pos = @{ $json->{$level}->{"ids"}->{$id}->{nonRegulator} };
            for my $i ( 0 .. $#Pos ) {
                my ( $absstart, $absend ) = @{ $Pos[$i] }[ 1, 2 ];
                my ( $start,    $end )    = @{ ${ $json->{$level}->{"ids"}->{$id}->{nonRegulatorCor} }[$i] };
                if ( $pos1 >= $absstart && $pos1 <= $absend ) {
                    my $relative = $start + ( $pos1 - $absstart ) * ( $end - $start ) / abs( $end - $start );
                    return $relative;
                }
            }
        }
        else {
            print STDERR "$id\t$level records is wrong\n";
            return "-1";
        }
    }
}

sub abs2rel {
    $json = shift @_;
    my $option_abs2rel = shift @_;
    chomp($option_abs2rel);
    while ( $option_abs2rel =~ /([^;]*)=([^;]*)/g ) {
        my ( $level, $cortype ) = ( $1, $2 );
        if ( $cortype eq "relative" ) {

            #level can only be product
            absolute2relative($level);
        }
        elsif ( $cortype eq "absolute" ) {

        }
        else {
            print STDERR "Unkonwn abs2rel type...exit\n";
            exit 255;
        }
    }
}

sub absolute2relative {
    my $level = shift @_;
    if ( exists $json->{$level} ) {
        my @ids_order = @{ $json->{$level}->{turns} };
        for my $id (@ids_order) {
            my ( $start, $end );
            if ( $json->{$level}->{"ids"}->{$id}->{"nonRegulatorCor"} ne "" ) {
                my @Loc = @{ $json->{$level}->{"ids"}->{$id}->{nonRegulatorCor} };
                for my $i ( 0 .. $#Loc ) {
                    my ( $start, $end ) = @{ $Loc[$i] };
                    my @values = @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"}[$i] };
                    if ( ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[$i]->[3] eq "+" ) {
                        ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[$i]->[1] = "-" . $start;
                        ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[$i]->[2] = "-" . $end;
                    }
                    else {
                        ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[$i]->[1] = "-" . $end;
                        ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[$i]->[2] = "-" . $start;
                    }
                }

            }
        }
    }
}

sub merger {
    $json = shift @_;
    our $sub_json = shift @_;

    my $ref_file        = shift @_;
    my $sub_file        = shift @_;
    my $overlap_percent = shift @_;
    my @header          = @_;
    if ( exists $sub_json->{locus}->{ids} ) {
        foreach my $id ( sort keys %{ $sub_json->{locus}->{ids} } ) {
            my $sub_chr    = $sub_json->{locus}->{ids}->{$id}->{seqId};
            my $sub_start  = $sub_json->{locus}->{ids}->{$id}->{min_start};
            my $sub_end    = $sub_json->{locus}->{ids}->{$id}->{max_end};
            my $sub_strand = $sub_json->{locus}->{ids}->{$id}->{strand};
            my $match_id;
            foreach my $ref_id ( sort keys %{ $json->{locus}->{ids} } ) {
                my $ref_chr    = $json->{locus}->{ids}->{$ref_id}->{seqId};
                my $ref_start  = $json->{locus}->{ids}->{$ref_id}->{min_start};
                my $ref_end    = $json->{locus}->{ids}->{$ref_id}->{max_end};
                my $ref_strand = $json->{locus}->{ids}->{$ref_id}->{strand};
                if ( $sub_chr eq $ref_chr && $sub_strand eq $ref_strand ) {
                    if ( $sub_start == $ref_start && $sub_end == $ref_end ) {
                        $match_id = $ref_id;
                        last;
                    }
                    else {
                        my %a;
                        $a{$sub_start} = 1;
                        $a{$sub_end}   = 1;
                        $a{$ref_start} = 2;
                        $a{$ref_end}   = 2;
                        my @c = sort { $a <=> $b } keys %a;
                        my @d = map  { $a{$_} } @c;
                        my $str     = join( "", @d );
                        my $overlap = $c[2] - $c[1];

                        if (   $str != 1122
                            && $str != 2211
                            && ( $overlap / ( $ref_end - $ref_start + 1 ) >= $overlap_percent || $overlap / ( $sub_end - $sub_start + 1 ) >= $overlap_percent )
                          )
                        {    #overlap is based on ref or sub
                            $match_id = $ref_id;
                            last;
                        }
                    }
                    ## -------------- ref    ref_end-sub_start
                    ##       --------------
                    #
                    ## --------------
                    ##       --------------- ref	sub_end-ref_start
                }
            }
            if ( defined $match_id ) {    ## have the same position
                print STDERR "ref $ref_file and sub $sub_file have close position locus ref:$match_id and sub:$id\n";
                restore_subjson( $id, "same", $ref_file, $sub_file, $match_id );
                renew_max_minPos( "locus", $match_id, $sub_start, $sub_end );

            }
            else {
                push @{ $json->{locus}->{turns} }, $id;
                $json->{locus}->{ids}->{$id} = $sub_json->{locus}->{ids}->{$id};
                print STDERR "ref $ref_file and sub $sub_file: sub has new locus id $id\n";
                restore_subjson( $id, "new", $ref_file, $sub_file );
            }

        }
    }
    if ( $#header != -1 ) {
        @{ $json->{"header"} } = @header;
    }

    sorter();

    sub renew_max_minPos {
        my ( $level, $id, $start, $end ) = @_;
        if ( !defined $json->{$level}->{"ids"}->{$id}->{"min_start"} ) {
            $json->{$level}->{"ids"}->{$id}->{"min_start"} = $start;
            $json->{$level}->{"ids"}->{$id}->{"max_end"}   = $end;
        }
        elsif ( defined $json->{$level}->{"ids"}->{$id}->{"min_start"} ) {
            if ( $json->{$level}->{"ids"}->{$id}->{"min_start"} > $start ) {
                $json->{$level}->{"ids"}->{$id}->{"min_start"} = $start;
            }
            if ( $json->{$level}->{"ids"}->{$id}->{"max_end"} < $end ) {
                $json->{$level}->{"ids"}->{$id}->{"max_end"} = $end;
            }
        }
        if ( $level eq "locus" ) {
            ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[0]->[1] = $json->{$level}->{"ids"}->{$id}->{"min_start"};
            ${ $json->{$level}->{"ids"}->{$id}->{nonRegulator} }[0]->[2] = $json->{$level}->{"ids"}->{$id}->{"max_end"};
        }
    }

    sub restore_subjson {
        my $id          = shift @_;
        my $type        = shift @_;
        my $ref_file    = shift @_;
        my $sub_file    = shift @_;
        my $ref_locusid = shift @_;
        my ( @ref_primaryids, @ref_processedids, @ref_productids );
        if ( defined $ref_locusid && exists $json->{locus}->{ids}->{$ref_locusid}->{childId} ) {
            @ref_primaryids = sort keys %{ $json->{locus}->{ids}->{$ref_locusid}->{childId} };
            for my $ref_primaryid (@ref_primaryids) {
                if ( exists $json->{primary}->{ids}->{$ref_primaryid}->{childId} ) {
                    push @ref_processedids, sort keys %{ $json->{primary}->{ids}->{$ref_primaryid}->{childId} };
                    for my $ref_processedid (@ref_processedids) {
                        if ( exists $json->{processed}->{ids}->{$ref_processedid}->{childId} ) {
                            push @ref_productids, sort keys %{ $json->{processed}->{ids}->{$ref_processedid}->{childId} };    ##debug:previous is product
                        }
                    }
                }
            }
        }
        if ( exists $sub_json->{locus}->{ids}->{$id}->{childId} ) {
            foreach my $childId_primary ( sort keys %{ $sub_json->{locus}->{ids}->{$id}->{childId} } ) {
                evaluate_samePos( $childId_primary, "primary", $type, $ref_file, $sub_file, $ref_locusid, @ref_primaryids );
                if ( exists $sub_json->{primary}->{ids}->{$childId_primary}->{childId} ) {
                    foreach my $childId_processed ( sort keys %{ $sub_json->{primary}->{ids}->{$childId_primary}->{childId} } ) {
                        evaluate_samePos( $childId_processed, "processed", $type, $ref_file, $sub_file, "", @ref_processedids );
                        if ( exists $sub_json->{processed}->{ids}->{$childId_processed}->{childId} ) {
                            foreach my $childId_product ( sort keys %{ $sub_json->{processed}->{ids}->{$childId_processed}->{childId} } ) {
                                evaluate_samePos( $childId_product, "product", $type, $ref_file, $sub_file, "", @ref_productids );
                            }
                        }
                    }
                }
            }
        }

    }

    sub evaluate_samePos {
        my ( $subid, $level, $type, $ref_file, $sub_file, $ref_locusid, @refids ) = @_;
        if ( $type eq "new" ) {
            push @{ $json->{$level}->{turns} }, $subid;
            $json->{$level}->{ids}->{$subid} = $sub_json->{$level}->{ids}->{$subid};
            return;
        }
        my $same = 0;
        for my $refid (@refids) {
            my $times = 0;
            if ( $#{ $json->{$level}->{"ids"}->{$refid}->{nonRegulator} } == $#{ $sub_json->{$level}->{"ids"}->{$subid}->{nonRegulator} } ) {
                foreach my $ind ( 0 .. $#{ $json->{$level}->{"ids"}->{$refid}->{nonRegulator} } ) {
                    if ( ${ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} }[$ind]->[0] eq
                           ${ $sub_json->{$level}->{"ids"}->{$subid}->{"nonRegulator"} }[$ind]->[0]
                        && ${ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} }[$ind]->[3] eq
                        ${ $sub_json->{$level}->{"ids"}->{$subid}->{"nonRegulator"} }[$ind]->[3]
                        && ${ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} }[$ind]->[1] ==
                        ${ $sub_json->{$level}->{"ids"}->{$subid}->{"nonRegulator"} }[$ind]->[1]
                        && ${ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} }[$ind]->[2] ==
                        ${ $sub_json->{$level}->{"ids"}->{$subid}->{"nonRegulator"} }[$ind]->[2] )
                    {
                        $times++;
                    }
                }
            }
            my $num = $#{ $json->{$level}->{"ids"}->{$refid}->{nonRegulator} } + 1;

            if ( $times == $num ) {
                print STDERR "inevaluate_samePos ref $ref_file and sub $sub_file have same position $level ref:$refid and sub:$subid. skip store $subid\n";
                $same = 1;

                my @all_values = @{ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} };
                delete $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"};
                foreach my $each (@all_values) {
                    my ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @{$each};
                    if ( $attributes =~ /(.*Dbxref=)(.*)/ ) {
                        $attributes = $1 . "$sub_file:$subid," . $2;
                    }
                    else {
                        $attributes = $attributes . "Dbxref=$sub_file:$subid;";
                    }
                    my @values = ( $seqId, $start, $end, $strand, $featureType, $attributes );
                    push @{ $json->{$level}->{"ids"}->{$refid}->{"nonRegulator"} }, \@values;
                }
                return;
            }
        }
        if ( !$same ) {
            push @{ $json->{$level}->{turns} }, $subid;
            $json->{$level}->{ids}->{$subid} = $sub_json->{$level}->{ids}->{$subid};

            if ( $level eq "primary" ) {
                $json->{locus}->{"ids"}->{$ref_locusid}->{"childId"}->{$subid} = 1000;
                delete $json->{$level}->{ids}->{$subid}->{parentId};
                $json->{$level}->{ids}->{$subid}->{parentId}->{$ref_locusid} = 1;
            }
            return;
        }
    }

}

sub detector_informat {
    my $infh = shift @_;
    my $informat;
    while (<$infh>) {
        chomp;
        next if /^#/ || /^\s*$/;
        my @line = split /\t/;
        if ( $line[8] =~ /gene_id\s+".*";\s+transcript_id\s+".*"/ ) {
            $informat = "gtf";
        }
        elsif ( $line[3] eq "+" || $line[3] eq "-" ) {
            $informat = "ngff";
        }
        elsif ( $line[8] =~ /=/ ) {
            $informat = "gff3";
        }
    }
    return $informat;
}

sub detector_outformat {
    my $outfile = shift @_;
    ( my $outformat = $outfile ) =~ s/(.*)\.(.*?)$/$2/;
    if ( $outformat ne "gff3" && $outformat ne "ngff" && $outformat ne "json" && $outformat ne "gtf" && $outformat ne "ntsv" ) {
        print STDERR "undefined output format($outformat)...exit\n";
        exit 255;
    }
    return $outformat;
}

sub isPossibleCombination {
    my ( $informat, $outformat ) = @_;
    my $Trans = {
        "ngff" => { "ngff" => 1, "gtf"  => 1, "gff3" => 1, "json" => 1, "ntsv" => 1 },
        "gtf"  => { "ngff" => 1, "gff3" => 1, "gtf"  => 1 },
        "gff3" => { "ngff" => 1, "gff3" => 1, "gtf"  => 1 },
    };
    if ( exists $Trans->{$informat}->{$outformat} ) {
        print STDERR "#Cheer, Input format($informat) can produce output format($outformat)\n";
    }
    else {
        print STDERR Dumper $Trans;
        print STDERR "#Wrong, Input format($informat) could not produce output format($outformat)\n";
        exit 255;
    }
}
##parsers
our $parsers = {
    ngff => \&parser_ngff,
    gff3 => \&parser_gff3,
    gtf  => \&parser_gtf,
    json => \&parser_json
};

sub parser_ngff {
    my $infh     = shift @_;
    my $prefix   = shift @_;
    my $isSimple = shift @_;
    $json = {};
    $json->{"input_format"} = "ngff";
    push @{ $json->{"header"} }, "";
    while (<$infh>) {
        chomp;
        next if /^\s*$/;
        if (/^#/) {
            push @{ $json->{"header"} }, $_;
        }
        else {
            my ( $seqId, $start, $end, $strand, $level, $featureType, $id, $parentId, $attributes ) = split /\t/;
            $id = $prefix . $id;
            if ( $level ne "chromosome" ) {
                $parentId =~ s/^/$prefix/g;
                $parentId =~ s/,/,$prefix/g;
            }
            my ( $attr_hash, $attr_tag ) = parse_Attr( $attributes, $json->{"input_format"}, $prefix );
            store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );
        }
    }
    sort_inSubStructure();
    store_childID();

    sorter();

    if ( !defined $isSimple ) {
        store_gid($prefix);    #store relation of gid and pid; store gene position(min_s max_e)
    }
    return $json;
}

sub store_gid {
    my $prefix = shift @_;
    my %hash;
    for my $id ( @{ $json->{product}->{"turns"} } ) {
        my $gene;
        if ( exists $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene} ) {
            $gene = $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene};

        }
        elsif ( exists $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene_id} ) {
            $gene = $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene_id};
        }
        else {
            print STDERR "Check sub(store_gid),exit...\n";
            exit 255;
        }

        if ( !exists $hash{$gene}{$id} ) {
            push @{ $Gene{gid2productid}{$gene}{childId} }, $id;
            $hash{$gene}{$id} = 1;
        }
        if ( exists $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene_name} ) {
            $Gene{gid2name}{$gene} = $json->{product}->{ids}->{$id}->{nonRegulator_shared_attributes}->{gene_name};
        }

    }
    undef %hash;
    foreach my $gene ( sort keys %{ $Gene{gid2productid} } ) {
        for my $id ( @{ $Gene{gid2productid}{$gene}{childId} } ) {
            foreach my $parentId (
                sort { $json->{processed}->{ids}->{$a}->{min_start} <=> $json->{processed}->{ids}->{$b}->{min_start} }
                keys %{ $json->{product}->{ids}->{$id}->{parentId} }
              )
            {    ## parentId is processed id
                if ( !exists $hash{$gene}{$id} ) {
                    push @{ $Gene{gid2processid}{$gene}{childId} }, $parentId;
                    $hash{$gene}{$id} = 1;
                }

                if ( !defined $Gene{gid2productid}{$gene}{min_start} ) {
                    $Gene{gid2productid}{$gene}{min_start} = $json->{processed}->{ids}->{$parentId}->{min_start};
                    $Gene{gid2productid}{$gene}{max_end}   = $json->{processed}->{ids}->{$parentId}->{max_end};
                }
                elsif ( defined $Gene{gid2productid}{$gene}{min_start} ) {
                    if ( $Gene{gid2productid}{$gene}{min_start} > $json->{processed}->{ids}->{$parentId}->{min_start} ) {
                        $Gene{gid2productid}{$gene}{min_start} = $json->{processed}->{ids}->{$parentId}->{min_start};
                    }
                    if ( $Gene{gid2productid}{$gene}{max_end} < $json->{processed}->{ids}->{$parentId}->{max_end} ) {
                        $Gene{gid2productid}{$gene}{max_end} = $json->{processed}->{ids}->{$parentId}->{max_end};
                    }
                }
            }
        }
    }
    foreach my $gene ( keys %{ $Gene{gid2processid} } ) {
        @{ $Gene{gid2processid}{$gene}{childId} } = &Uniq_array( @{ $Gene{gid2processid}{$gene}{childId} } );
    }
}

sub store_childID {
    my @levels_order = ( "chromosome", "locus", "primary", "processed", "product" );
    for my $i ( 0 .. $#levels_order ) {
        my $level = $levels_order[$i];
        if ( exists $json->{$level} ) {
            my $num = 1;
            for my $id ( @{ $json->{$level}->{"turns"} } ) {
                Coordinate( $id, $level );
                foreach my $parentId ( keys %{ $json->{$level}->{"ids"}->{$id}->{parentId} } ) {
                    next if ( !exists $Levels->{$parentId} );
                    my $parentLevel = $Levels->{$parentId};
                    $json->{$parentLevel}->{"ids"}->{$parentId}->{"childId"}->{$id} = $num;
		    inheritance_implementation($id,$parentId,$level,$parentLevel);
                }
                $num++;
            }
        }
    }
}

my $inheritable = {'edited_site'=>1};
sub inheritance_implementation{
  my ($id,$parentId,$level,$parentLevel)=@_;
  my $nonRegulator_shared_attributes_hash=$json->{$parentLevel}->{"ids"}->{$parentId}->{"nonRegulator_shared_attributes"};
  foreach my $key (keys %{$nonRegulator_shared_attributes_hash}){
	if(exists $inheritable->{$key}){
		my $value=$nonRegulator_shared_attributes_hash->{$key};
		$json->{$level}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{$key}=$value;
	}
  }
}


sub store_common {
    my ( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId ) = @_;
    if ( !exists $Levels->{$id} ) {
        push @{ $json->{$level}->{"turns"} }, $id;
        $Levels->{$id} = $level;
    }

    $json->{$level}->{"ids"}->{$id}->{"seqId"} = $seqId;
    if ( !defined $json->{$level}->{"ids"}->{$id}->{"min_start"} && $start > 0 && $featureType ne "regulator" ) {
        $json->{$level}->{"ids"}->{$id}->{"min_start"} = $start;
        $json->{$level}->{"ids"}->{$id}->{"max_end"}   = $end;
    }
    elsif ( defined $json->{$level}->{"ids"}->{$id}->{"min_start"} && $start > 0 && $featureType ne "regulator" ) {
        if ( $json->{$level}->{"ids"}->{$id}->{"min_start"} < $start ) {
            $json->{$level}->{"ids"}->{$id}->{"min_start"} = $json->{$level}->{"ids"}->{$id}->{"min_start"};
        }
        else {
            $json->{$level}->{"ids"}->{$id}->{"min_start"} = $start;
        }
        if ( $json->{$level}->{"ids"}->{$id}->{"max_end"} > $end ) {
            $json->{$level}->{"ids"}->{$id}->{"max_end"} = $json->{$level}->{"ids"}->{$id}->{"max_end"};
        }
        else {
            $json->{$level}->{"ids"}->{$id}->{"max_end"} = $end;
        }
    }
    elsif ( $start < 0 && $end < 0 && $level eq "product" ) {
        $json->{$level}->{"ids"}->{$id}->{"min_start"} = $start;
        $json->{$level}->{"ids"}->{$id}->{"max_end"}   = $end;
    }
    $json->{$level}->{"ids"}->{$id}->{"strand"}      = $strand;
    $json->{$level}->{"ids"}->{$id}->{"featureType"} = $featureType;
    my $featureUid = "${seqId}_${start}_${end}_${strand}";

    delete $attr_hash->{ID};             #not shared attr
    delete $attr_hash->{exon_number};    #not shared attr
    delete $attr_hash->{exon_id};        #not shared attr

    if ( $attributes eq "." ) {
        $attributes = "featureUid=$featureUid;";
    }
    elsif ( !exists $attr_hash->{featureUid} ) {
        $attributes = $attributes . "featureUid=$featureUid;";
    }
    my @values = ( $seqId, $start, $end, $strand, $featureType, $attributes );

    if ( $featureType eq "regulator" || $featureType eq "polyA_sequence" ) {
        my $regulator_biotype;
        if ( exists $attr_hash->{regulator_biotype} ) {
            $regulator_biotype = $attr_hash->{regulator_biotype};
        }
        else {
            $regulator_biotype = $featureType;
        }
        push @{ $json->{$level}->{"ids"}->{$id}->{"regulator"}->{$regulator_biotype} }, \@values;
        $json->{$level}->{"ids"}->{$id}->{"regulator_shared_attributes"} = $attr_hash;
        my $exists_shared_attributes=$json->{$level}->{"ids"}->{$id}->{"regulator_shared_attributes"};
	my $new_shared_attributes=add_shared_attributes($exists_shared_attributes,$attr_hash);
        $json->{$level}->{"ids"}->{$id}->{"regulator_shared_attributes"} = $new_shared_attributes;
    }
    else {
        push @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }, \@values;
        my $exists_shared_attributes=$json->{$level}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"};
	my $new_shared_attributes=add_shared_attributes($exists_shared_attributes,$attr_hash);
        $json->{$level}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"} = $new_shared_attributes;
    }
sub add_shared_attributes{
my ($exists_shared_attributes,$attr_hash)=@_;
my $new_shared_attributes;
	if (!defined $exists_shared_attributes){
		$new_shared_attributes=$attr_hash;
	}else{
		foreach my $key (keys %{$attr_hash}){
			my $new_value=$attr_hash->{$key};
			if (exists $exists_shared_attributes->{$key} && $exists_shared_attributes->{$key} ne $attr_hash->{$key} ){
				$exists_shared_attributes->{$key}=$exists_shared_attributes->{$key}.",$new_value";
			}elsif(!exists $exists_shared_attributes->{$key}){
				$exists_shared_attributes->{$key}=$new_value;
			}
		}
		$new_shared_attributes=$exists_shared_attributes;
	}
}

    if ( defined $parentId ) {
        my @parentIds = split /,/, $parentId;
        map { $json->{$level}->{"ids"}->{$id}->{parentId}->{$_} = 1 } @parentIds;
    }
    if ( $level eq "processed" && ( exists $json->{$level}->{"ids"}->{$id}->{nonRegulator_shared_attributes}->{gene_id} ) && $json->{input_format} ne "gff3" ) {
        my $gid = $json->{$level}->{"ids"}->{$id}->{nonRegulator_shared_attributes}->{gene_id};

        $Gene{tid2gid}{$id} = $gid;
        push @{ $Gene{gid2processid}{$gid}{childId} }, $id;
        $Gene{gidstart}{$gid} = $json->{locus}->{"ids"}->{ $gid . ".locus" }->{"min_start"};
        $Gene{gidend}{$gid}   = $json->{locus}->{"ids"}->{ $gid . ".locus" }->{"max_end"};
    }
    elsif ( $level eq "processed" && ( exists $json->{$level}->{"ids"}->{$id}->{nonRegulator_shared_attributes}->{gene} ) && $json->{input_format} ne "gff3" ) {
        my $gid = $json->{$level}->{"ids"}->{$id}->{nonRegulator_shared_attributes}->{gene};
        $Gene{tid2gid}{$id} = $gid;
        push @{ $Gene{gid2processid}{$gid}{childId} }, $id;
        my $start = $json->{"processed"}->{"ids"}->{$id}->{"min_start"};
        my $end   = $json->{"processed"}->{"ids"}->{$id}->{"max_end"};
        if ( !exists $Gene{gidstart}{$gid} || $start < $Gene{gidstart}{$gid} ) {
            $Gene{gidstart}{$gid} = $start;
        }
        if ( !exists $Gene{gidend}{$gid} || $end > $Gene{gidend}{$gid} ) {
            $Gene{gidend}{$gid} = $end;
        }
    }
    elsif ( $level eq "processed" && $json->{input_format} eq "gff3" ) {
        my $pri_id = ( keys %{ $json->{$level}->{ids}->{$id}->{parentId} } )[0];
        my $gid    = $json->{primary}->{ids}->{$pri_id}->{nonRegulator_shared_attributes}->{Parent};

        $Gene{tid2gid}{$id} = $gid;
        push @{ $Gene{gid2processid}{$gid}{childId} }, $id;
        my $start = $json->{"processed"}->{"ids"}->{$id}->{"min_start"};
        my $end   = $json->{"processed"}->{"ids"}->{$id}->{"max_end"};
        if ( !exists $Gene{gidstart}{$gid} || $start < $Gene{gidstart}{$gid} ) {
            $Gene{gidstart}{$gid} = $start;
        }
        if ( !exists $Gene{gidend}{$gid} || $end > $Gene{gidend}{$gid} ) {
            $Gene{gidend}{$gid} = $end;
        }
    }
    elsif ( $level eq "processed" && $json->{input_format} ne "gff3" ) {
        my $gid = ( keys %{ $json->{$level}->{ids}->{$id}->{parentId} } )[0];
        $Gene{gidstart}{$gid} = $json->{"primary"}->{"ids"}->{$gid}->{"min_start"};
        $Gene{gidend}{$gid}   = $json->{"primary"}->{"ids"}->{$gid}->{"max_end"};
    }
}

sub CDSphase {
    foreach my $id ( keys %{ $json->{product}->{"ids"} } ) {
        my %hash_tmp;
        my @Values    = @{ $json->{product}->{"ids"}->{$id}->{"nonRegulator"} };
        my @ValuesCor = @{ $json->{product}->{"ids"}->{$id}->{"nonRegulatorCor"} };
        foreach my $i ( 0 .. $#ValuesCor ) {
            my $each = $ValuesCor[$i];
            my ( $start, $end ) = @$each;
            $hash_tmp{$start}{end}  = $end;
            $hash_tmp{$start}{indx} = $i;
        }
        my $last_mod;
        my $phase;
        my $num = 0;
        foreach my $rstart ( sort { $a <=> $b } keys %hash_tmp ) {
            my $indx = $hash_tmp{$rstart}{indx};
            my $each = $Values[$indx];
            my ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$each;
            next if ( $featureType ne "CDS" || $start < 0 );
            if ( $num == 0 ) {
                $last_mod = ( $end - $start + 1 ) % 3;
                $phase    = 0;
            }
            else {
                $phase = ( 3 - $last_mod ) % 3;
                my $frameshift = 0;
                if ( $attributes =~ /frameshift=(\d+):-(\d+)/ ) {
                    $frameshift = $2;
                }
                $last_mod = ( ( $end - $start + 1 ) + $frameshift - $phase ) % 3;
            }

            $Gene{CDSphase}{$id}{$indx} = $phase;
            if ( $attributes =~ /phase=(\d)/ ) {
                my $phase_old = $1;
                if ( $phase_old ne $phase ) {
                    print STDERR "$id($start - $end)'s phase is wrong, been corrected from $phase_old to $phase\n";
                }
                $attributes =~ s/phase=.;//;
            }
            $each->[5] = $attributes;
            $num++;
        }

        $num = 0;
        for my $each ( @{ $json->{product}->{"ids"}->{$id}->{"nonRegulator"} } ) {
            my ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$each;
            next if ( $featureType ne "CDS" || $start < 0 );
            my $phase = $Gene{CDSphase}{$id}{$num};
            if ( !defined $phase ) {
                print STDERR "undefined phase in attributes:$attributes\n";
            }
            $each->[5] = $attributes . "phase=$phase;";
            $num++;

        }

    }
}

sub Coordinate {
    my ( $last_id, $last_level ) = @_;

    my @S_sub;
    my @E_sub;
    my @Strand_sub;
    my %fusion_hash;

    if ( exists $json->{$last_level}->{"ids"}->{$last_id}->{"nonRegulator"} ) {
        foreach my $i ( 0 .. $#{ $json->{$last_level}->{"ids"}->{$last_id}->{"nonRegulator"} } ) {
            my $each = ${ $json->{$last_level}->{"ids"}->{$last_id}->{"nonRegulator"} }[$i];
            my ( $seqIdl, $startl, $endl, $strandl, $featureType, $attributesl, $fusion_type ) = @$each;
            if ( !defined $fusion_type ) {
                $fusion_type = 1;
            }
            if ( $strandl eq "-" && $startl > $endl ) {    #when strand is -, start also need to be less than end
                my $tmp = $startl;
                $startl = $endl;
                $endl   = $tmp;
            }
            push @{ $fusion_hash{$fusion_type}{S} },      $startl;
            push @{ $fusion_hash{$fusion_type}{E} },      $endl;
            push @{ $fusion_hash{$fusion_type}{Strand} }, $strandl;
        }
    }

    foreach my $fusion_type ( sort { $a <=> $b } keys %fusion_hash ) {
        my @Strand_sub = @{ $fusion_hash{$fusion_type}{Strand} };
        if ( $Strand_sub[0] eq "+" ) {
            my @S_sub = @{ $fusion_hash{$fusion_type}{S} };
            my @E_sub = @{ $fusion_hash{$fusion_type}{E} };
            @{ $fusion_hash{$fusion_type}{S} }      = reverse(@S_sub);
            @{ $fusion_hash{$fusion_type}{E} }      = reverse(@E_sub);
            @{ $fusion_hash{$fusion_type}{Strand} } = reverse(@Strand_sub);
        }
    }

    # for relative abs2rels
    my $length_nonRegulator = 0;
    my @relLoc;
    foreach my $fusion_type ( sort { $a <=> $b } keys %fusion_hash ) {
	my @relLoc_sub;
        my @Strand_sub = @{ $fusion_hash{$fusion_type}{Strand} };
        my @S_sub = @{ $fusion_hash{$fusion_type}{S} };
        my @E_sub = @{ $fusion_hash{$fusion_type}{E} };
        for ( my $i = $#S_sub ; $i >= 0 ; $i = $i - 1 ) {
            my ( $this_start, $this_end );
            if ( $Strand_sub[$i] eq "+" ) {
                $this_start = $length_nonRegulator + 1;
                $this_end   = $this_start + abs( $E_sub[$i] - $S_sub[$i] );
            }
            else {
                $this_end   = $length_nonRegulator + 1;
                $this_start = $this_end + abs( $E_sub[$i] - $S_sub[$i] );
            }
            my @loc = ( $this_start, $this_end );
            if ( $Strand_sub[0] eq "+" ) {
                push @relLoc_sub, \@loc;
            }
            else {
                unshift @relLoc_sub, \@loc;
            }
            $length_nonRegulator += abs( $E_sub[$i] - $S_sub[$i] + 1 );
        }
	@relLoc=(@relLoc,@relLoc_sub);
    }
    @{ $json->{$last_level}->{"ids"}->{$last_id}->{"nonRegulatorCor"} } = @relLoc;

    my ( $pos, $abs_start, $rel_start, $rel_end ) = @_;
    if ( exists $json->{$last_level}->{"ids"}->{$last_id}->{"regulator"} ) {
        foreach my $regulator_biotype ( keys %{ $json->{$last_level}->{"ids"}->{$last_id}->{"regulator"} } ) {
            for my $each ( @{ $json->{$last_level}->{"ids"}->{$last_id}->{"regulator"}->{$regulator_biotype} } ) {
                my ( $seqIdl, $startl, $endl, $strandl, $featureType, $attributesl ) = @$each;
                my @loc = ( $startl, $endl );
                push @{ $json->{$last_level}->{"ids"}->{$last_id}->{"regulatorCor"}->{$regulator_biotype} }, \@loc;
            }
        }
    }

}

sub parse_Attr {
    my $str       = shift @_;
    my $InputType = shift @_;
    my $prefix    = shift @_;
    my $level     = shift @_;
    my @turns     = ();
    my %hash;
    if ( $InputType eq "gff3" || $InputType eq "ngff" || $InputType eq "gtf" ) {
        my @tmp;
        if ( $InputType eq "gtf" ) {
            $str =~ s/"//g;
            $str =~ s/;//g;
            @tmp = split /\s+/, $str;    ### bug: In this way, the value could not have space.
        }
        else {
            @tmp = split /;|=/, $str;
        }
        for ( my $i = 0 ; $i < $#tmp ; $i += 2 ) {
            my $tag   = $tmp[$i];
            my $value = $tmp[ $i + 1 ];
            $hash{$tag} = $value;
            if ( $tag ne "ID" && $tag ne "Parent" && $tag ne "" ) {
                push @turns, $tag;
            }
        }
        if ( exists $hash{ID} ) {
            $hash{ID} = $prefix . $hash{ID};
        }
        if ( exists $hash{Parent} ) {
            $hash{Parent} = $prefix . $hash{Parent};
        }
        if ( exists $hash{protein_id} ) {
            $hash{protein_id} = $prefix . $hash{protein_id};
        }

        if ( $InputType eq "ngff" && exists $hash{gene} ) {
            $hash{gene} = $prefix . $hash{gene};
        }
        if ( exists $hash{gene_id} ) {
            $hash{gene_id} = $prefix . $hash{gene_id};
        }
        if ( exists $hash{transcript_id} ) {
            $hash{transcript_id} = $prefix . $hash{transcript_id};
        }
    }
    elsif ( $InputType eq "gene" ) {
        while ( $str =~ /(\S+?)=(\S+?);/g ) {
            my $tag   = $1;    # why need redefined?
            my $value = $2;
            if ( $tag =~ /gene/ ) {
                $hash{$tag} = $value;
                push @turns, $tag;
            }
        }
    }
    return ( \%hash, \@turns );
}

sub AttrHash2String {
    my $attr_hash   = shift @_;
    my $attr_tag    = shift @_;
    my $featureType = shift @_;
    my $id          = shift @_;
    my $parentId    = shift @_;
    my $attributes  = "";
    my @Add;
    if ( $featureType ne "reference" ) {
        @Add = @{$attr_tag};
    }
    if ( $featureType eq "reference" && !exists $attr_hash->{seqSource} ) {
        $attr_hash->{seqSource} = "unknown";
        push @Add, "seqSource";
    }
    if ( $featureType eq "reference" && !exists $attr_hash->{circular} ) {
        $attr_hash->{circular} = "ND";
        push @Add, "circular";
    }
    if ( $featureType eq "region" && !exists $attr_hash->{primary_biotype} ) {
        $attr_hash->{primary_biotype} = "primary_transcript";
        push @Add, "primary_biotype";
    }
    if ( $featureType eq "region" && !exists $attr_hash->{transcribed} ) {
        $attr_hash->{transcribed} = "TRUE";
        push @Add, "transcribed";
    }
    if ( $featureType eq "cap" || $featureType eq "exon" || $featureType eq "intron" || $featureType eq "polyA_sequence" ) {
        if ( !exists $attr_hash->{processed_biotype} ) {
            $attr_hash->{processed_biotype} = "processed_transcript";
            push @Add, "processed_biotype";
        }
        if ( !exists $attr_hash->{product_biotype} && defined $parentId ) {
            my $gid = $json->{primary}->{ids}->{$parentId}->{nonRegulator_shared_attributes}->{Parent};
            if ( exists $Gene{type}{$gid} ) {
                $attr_hash->{product_biotype} = $Gene{type}{$gid};
            }
            else {
                $attr_hash->{product_biotype} = "unknown";

            }
            push @Add, "product_biotype";
        }

    }
    if ( $featureType eq "CDS" || $featureType eq "mature" ) {
        if ( !exists $attr_hash->{product_biotype} ) {
            if ( defined $parentId ) {
                my $product_biotype = $json->{processed}->{ids}->{$parentId}->{nonRegulator_shared_attributes}->{product_biotype};
                $attr_hash->{product_biotype} = $product_biotype;
            }
            else {
                $attr_hash->{product_biotype} = "unknown";
            }
        }
        push @Add, "product_biotype";

    }
    my %tmphash;
    @Add = grep { ++$tmphash{$_} < 2 } @Add;
    for my $each (@Add) {
        if ( exists $attr_hash->{$each} && defined $attr_hash->{$each} ) {
            if ( $attr_hash->{$each} =~ / / ) {
                $attributes .= "$each=\"$attr_hash->{$each}\";";
            }
            else {
                $attributes .= "$each=$attr_hash->{$each};";
            }
        }
    }
    return $attributes;
}

sub AddAttr {
    my ( $attr_hash, $attr_tag, $key, $value ) = @_;
    push @{$attr_tag}, $key;
    $attr_hash->{$key} = $value;
}

sub parser_gtf {
    my $infh   = shift @_;
    my $prefix = shift @_;
    $json = {};
    $json->{"input_format"} = "gtf";
    push @{ $json->{"header"} }, "";
    while (<$infh>) {
        chomp;
        next if /^\s*$/;
        if (/^#/) {
            push @{ $json->{"header"} }, $_;
        }
        else {
            my ( $seqId, undef, $featureType, $start, $end, undef, $strand, $phase, $attr ) = split /\t/;
            my ( $attr_hash, $attr_tag ) = parse_Attr( $attr, $json->{"input_format"}, $prefix );    # change
            my ( $parentId, $id, $protein_id, $product_id, $level, $gene_biotype, $gname, $gid, $tid );
            $gid = $attr_hash->{gene_id};
            $tid = $attr_hash->{transcript_id};
            if ( exists $attr_hash->{protein_id} ) {
                $protein_id = $attr_hash->{protein_id};
            }
            if ( exists $attr_hash->{gene_biotype} ) {
                $gene_biotype = $attr_hash->{gene_biotype};
            }
            else {
                $gene_biotype = "unknown";
            }
            if ( $featureType eq "transcript" ) {
                $level = "chromosome";
                my $attributes = AttrHash2String( $attr_hash, $attr_tag, "reference" );
                store_common( $seqId, $level, $seqId, $start, $end, "+", "reference", $attributes, $attr_hash, "genome" );
                $level = "locus";
                $id    = $gid . ".locus";
                if ( $gene_biotype eq "lincRNA" ) {
                    $gene_biotype = "lncRNA";
                }
                elsif ( $gene_biotype eq "protein_coding" ) {
                    $gene_biotype = "ORF";
                }
                elsif ( $gene_biotype eq "misc_RNA" ) {
                    $gene_biotype = "unknown";
                }
                $Gene{type}{$gid} = $gene_biotype;
                $attributes = AttrHash2String( $attr_hash, $attr_tag, "locus" );
                store_common( $id, $level, $seqId, $start, $end, $strand, "locus", $attributes, $attr_hash, $seqId );
                $level               = "primary";
                $id                  = $tid . ".primary";
                $Gene{tid2gid}{$tid} = $gid;
                $attributes          = AttrHash2String( $attr_hash, $attr_tag, "region" );
                store_common( $id, $level, $seqId, $start, $end, $strand, "region", $attributes, $attr_hash, $gid . ".locus" );
            }
            elsif ( $featureType eq "exon" ) {
                $id       = $tid;
                $parentId = $tid;
                $level    = "processed";
                my $gid        = $Gene{tid2gid}{$id};
                my $attributes = AttrHash2String( $attr_hash, $attr_tag, "exon" );
                store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId . ".primary" );

                if ( $Gene{type}{$gid} eq "lncRNA" ) {
                    $product_id  = $id . ".mature";
                    $level       = "product";
                    $featureType = "mature";

                    if ( !exists $attr_hash->{product_biotype} ) {
                        AddAttr( $attr_hash, $attr_tag, "product_biotype", "lncRNA" );
                    }
                    my $attributes = AttrHash2String( $attr_hash, $attr_tag, "mature" );
                    store_common( $product_id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );
                }
            }
            elsif ( $featureType eq "CDS" ) {
                $id       = $protein_id;
                $parentId = $tid;
                $level    = "product";
                my $gid = $Gene{tid2gid}{$parentId};
                if ( !exists $attr_hash->{product_biotype} ) {
                    AddAttr( $attr_hash, $attr_tag, "product_biotype", "lncRNA" );
                }
                if ( !exists $attr_hash->{gene} ) {
                    AddAttr( $attr_hash, $attr_tag, "gene", $gid );
                }
                my $attributes = AttrHash2String( $attr_hash, $attr_tag, "CDS" );
                store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );
            }
            else {
                next;
            }
        }
    }
    print STDERR "check: in gtf\n";
    sort_inSubStructure();
    redefined_Level("locus");
    store_childID();

    sorter();

    store_gid($prefix);
    return $json;
}

sub redefined_Level {
    my $level = shift @_;
    foreach my $id ( keys %{ $json->{chromosome}->{ids} } ) {    # no sort needed
        my $max_end = $json->{chromosome}->{ids}->{$id}->{max_end};
        my @Pos     = @{ $json->{chromosome}->{ids}->{$id}->{nonRegulator} };
        my $attr    = @{ $Pos[0] }[5];
        delete $json->{chromosome}->{ids}->{$id};
        $json->{chromosome}->{ids}->{$id}->{max_end}     = $max_end;
        $json->{chromosome}->{ids}->{$id}->{min_start}   = 1;
        $json->{chromosome}->{ids}->{$id}->{strand}      = "+";
        $json->{chromosome}->{ids}->{$id}->{featureType} = "reference";
        $json->{chromosome}->{ids}->{$id}->{seqId}       = $id;

        my @values = ( $id, "1", $max_end, "+", "reference", $attr );
        push @{ $json->{chromosome}->{ids}->{$id}->{nonRegulator} }, \@values;
        $json->{chromosome}->{ids}->{$id}->{parentId}->{genome} = 1;
    }
    if ( defined $level && $level eq "locus" ) {
        foreach my $id ( keys %{ $json->{locus}->{ids} } ) {    # not need sort
            my $end   = $json->{$level}->{ids}->{$id}->{max_end};
            my $start = $json->{$level}->{ids}->{$id}->{min_start};
            ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[0]->[1] = $start;
            ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[0]->[2] = $end;

            my ( $attr_hash, $attr_tag ) = parse_Attr( ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[0]->[5], "gene", "" );    # need change
            $json->{$level}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"} = $attr_hash;                                           # need change
            my $attr = "";
            foreach my $tag ( @{$attr_tag} ) {                                                                                          # need change
                $attr .= "$tag=$attr_hash->{$tag};";
            }
            ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[0]->[5] = $attr;
            my $final = ${ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }[0];
            delete $json->{$level}->{"ids"}->{$id}->{"nonRegulator"};
            push @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} }, $final;
        }
    }

}

sub parser_gff3 {
    my $infh    = shift @_;
    my $prefix  = shift @_;
    my $ncbi    = shift @_;
    my $haveChr = 0;
    $json = {};
    $json->{"input_format"} = "gff3";
    push @{ $json->{"header"} }, "";

    my $filepos = 0;

    while (<$infh>) {
        chomp;
        next if /^\s*$/;
        if (/^#/) {
            push @{ $json->{"header"} }, $_;
            $filepos = tell($infh);
        }
        else {
            last;
        }
    }
    seek( $infh, $filepos, 0 );
    while (<$infh>) {
        chomp;
        next if /^#/;
        my ( $seqId, undef, $featureType, $start, $end, undef, $strand, $phase, $attr ) = split /\t/;
        if ( $featureType eq "CDS" && $phase =~ /0-2/ ) {
            if ( $attr =~ /;$/ ) {
                $attr .= "phase=$phase;";
            }
            else {
                $attr .= ";phase=$phase;";

            }
        }
        my ( $attr_hash, $attr_tag ) = parse_Attr( $attr, $json->{"input_format"}, $prefix );    # change
        my ( $parentId, $id, $protein_id, $product_id, $level, $gene_biotype );
        my $gname = "unknown";
        $id = $attr_hash->{ID};
        if ( exists $attr_hash->{Parent} ) {
            $parentId = $attr_hash->{Parent};
        }
        elsif ( !$attr_hash->{Parent} && $featureType ne "chromosome" && $featureType ne "region" && $featureType ne "gene" && $featureType ne "pseudogene" ) {
            print STDERR "#WARNING: no parentId\t$id\n";
            next;
        }
        if ( exists $attr_hash->{protein_id} ) {
            $protein_id = $attr_hash->{protein_id};
        }
        if ( exists $attr_hash->{product} ) {
            $product_id = $attr_hash->{product};
        }
        if ( $featureType eq "chromosome" || $featureType eq "region" ) {
            $haveChr = 1;
            $level   = "chromosome";

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "reference" );
            store_common( $seqId, $level, $seqId, $start, $end, "+", "reference", $attributes, $attr_hash, "genome" );

        }
        elsif ( $featureType =~ "gene" || $featureType =~ "pseudogene" ) {
            my $gid = $id;

            if ( exists $attr_hash->{gene_biotype} ) {
                $gene_biotype = $attr_hash->{gene_biotype};
            }
            else {
                $gene_biotype = "unknown";
            }
            if ( exists $attr_hash->{Name} ) {
                $gname = $attr_hash->{Name};
                $attr_hash->{gene_name} = $gname;
            }
            else {
                $gname = "unknown";
            }
            if ( $gene_biotype eq "lincRNA" ) {
                $gene_biotype = "lncRNA";
            }
            elsif ( $gene_biotype eq "protein_coding" ) {
                $gene_biotype = "ORF";
            }
            elsif ( $gene_biotype eq "misc_RNA" ) {
                $gene_biotype = "unknown";
            }

            $Gene{type}{$gid}     = $gene_biotype;
            $Gene{gid2name}{$gid} = $gname;
            if ( !$haveChr ) {
                $level = "chromosome";

                my $attributes = AttrHash2String( $attr_hash, $attr_tag, "reference" );

                store_common( $seqId, $level, $seqId, $start, $end, "+", "reference", $attributes, $attr_hash, "genome" );

            }
            $level = "locus";
            my $id_Locus = $gid . ".locus";

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "locus" );
            store_common( $id_Locus, $level, $seqId, $start, $end, $strand, "locus", $attributes, $attr_hash, $seqId );

        }
        elsif (
               $featureType eq "mRNA"
            || $featureType eq "transcript"
            || (
                $featureType eq "primary_transcript"
                && (
                    (
                        exists $attr_hash->{product}
                        && ( $attr_hash->{product} =~ /microRNA/ || $attr_hash->{product} =~ /miRNA/ || $attr_hash->{product} =~ /miR\d+/ ) || $id =~ /miR\d+/
                    )
                )
            )

            || $featureType eq "lnc_RNA"
            || $featureType eq "tRNA" || $featureType eq "snRNA"
          )
        {
            $level = "primary";
            my $tid = $id;
            $id = $tid . ".primary";
            my $gid      = $parentId;
            my $id_Locus = $gid . ".locus";
            push @{ $Gene{gid2processid}{$gid}{childId} }, $tid;

            $Gene{tid2gid}{$tid} = $gid;

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "region" );
            store_common( $id, $level, $seqId, $start, $end, $strand, "region", $attributes, $attr_hash, $id_Locus );

        }
        elsif ( $featureType eq "exon" ) {
            $id    = $parentId;
            $level = "processed";
            my $gid;
            if ( exists $Gene{tid2gid}{$id} ) {
                $gid = $Gene{tid2gid}{$id};
                if ( !exists $Gene{type}{$gid} ) {
                    $Gene{type}{$gid} = "unknown";
                    print STDERR "#WARNING: $id\'s parent $gid has no biotype\n";
                }
            }
            else {
                print STDERR "#WARNING: ignore $id\n";
                next;
            }

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "exon", $id, $parentId . ".primary" );

            store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId . ".primary" );

            if ( $Gene{type}{$gid} eq "lncRNA" ) {
                $product_id  = $id . ".lncRNA.product";
                $level       = "product";
                $featureType = "mature";

                my $attributes = AttrHash2String( $attr_hash, $attr_tag, "mature", $product_id, $parentId );

                store_common( $product_id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );

            }
        }
        elsif ( $featureType eq "CDS" ) {
            my $gid;
            if ( defined $protein_id ) {
                $id = $protein_id;
            }
            else {
                ($id) = $parentId . ".product";
            }
            $level = "product";
            if ( exists $Gene{tid2gid}{$parentId} ) {
                $gid = $Gene{tid2gid}{$parentId};
                $attr_hash->{gene} = $gid;
            }
            else {
                print STDERR "#WARNING: ignore $parentId\n";
                next;
            }

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "CDS", $id, $parentId );
            store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );

        }
        elsif ( $featureType eq "miRNA" ) {

            my $gid = $Gene{tid2gid}{$parentId};

            if ( defined $ncbi && $ncbi eq "ncbi" && defined $product_id ) {
                ($id) = $product_id;
            }
            else {
                $id = $attr_hash->{ID} . ".miRNA.product";    # modified 0304
            }

            $level = "product";

            $featureType = "mature";

            my $attributes = AttrHash2String( $attr_hash, $attr_tag, "mature", $id );
            store_common( $id, $level, $seqId, $start, $end, $strand, $featureType, $attributes, $attr_hash, $parentId );

        }
        else {
            next;
        }
    }

    sort_inSubStructure();
    redefined_Level();
    store_childID();

    sorter();

    return $json;
}

sub sort_inSubStructure {
    my @Levels = ( "processed", "product" );
    foreach my $level (@Levels) {
	    next if !exists $json->{$level};
        my @ids = @{ $json->{$level}->{"turns"} };
        foreach my $id (@ids) {
            my %fusion_hash;

            my @Values = @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} };

            foreach my $i ( 0 .. $#Values ) {
                my $each = $Values[$i];
                my ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$each;
                my $fusion_type;
                if ( $attributes =~ /order=([^;]+)/ ) {
                    $fusion_type = $1;
                }
                else {
                    $fusion_type = 1;
                }
                push @{ $fusion_hash{$fusion_type} }, $i;
                $Values[$i]->[6] = $fusion_type;
            }

            my @Values_new;
            foreach my $fusion_type ( sort { $a <=> $b } keys %fusion_hash ) {
                my @Values_sub = @Values[ @{ $fusion_hash{$fusion_type} } ];
                if ( $Values_sub[0]->[1] > $Values_sub[-1]->[1] ) {
                    @Values_sub = reverse(@Values_sub);
                }
                @Values_new = ( @Values_new, @Values_sub );
            }

            @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} } = @Values_new;
        }
    }
}

sub parser_json { }

##autofillers
#our $autofillers = {
#    gff3 => \&autofiller_gff3,
#    gtf  => \&autofiller_gtf,
#    json => \&autofiller_json,
#    ntsv => \&autofiller_ntsv,
#    ngff => \&autofiller_ngff
#};
#
#sub autofiller_ngff {
#    $json = shift @_;
#
#}
#sub autofiller_gff3 { }
#sub autofiller_gtf  { }
#sub autofiller_json { }
#sub autofiller_ntsv { }

##printers
our $printers = {
    gff3 => \&printer_gff3,
    gtf  => \&printer_gtf,
    json => \&printer_json,
    ntsv => \&printer_ntsv,
    ngff => \&printer_ngff
};

sub printer_ngff {
    my ( $json, $outfh,$grep ) = @_;
    if (!defined $grep){
    CDSphase();
    }

    my %Done;
    my @Headers = @{ $json->{"header"} };
    my $flag    = 1;
    if ( $#Headers > -1 ) {
	print $outfh "##NGFF version 1.0\n";
        for my $each (@Headers) {
            next if $each eq "";
            print $outfh "$each\n";
            if ( $each =~ /^#seqId\tstart\tend\tstrand\tlevel\tfeatureType\tid\tparentId\tattributes/ ) {
                $flag = 0;
            }
        }
    }
    if ($flag) {
        print $outfh "#seqId\tstart\tend\tstrand\tlevel\tfeatureType\tid\tparentId\tattributes\n";
    }
    my ( $seqId, $start, $end, $strand, $level, $featureType, $id, $parentId, $attributes );
    foreach my $seqId ( @{ $json->{chromosome}->{turns} } ) {
        printer_subngff( $seqId, "chromosome", $outfh );
    }

    foreach my $locusId ( @{ $json->{locus}->{turns} } ) {
        printer_subngff( $locusId, "locus", $outfh );
        foreach my $priId (
            sort { $json->{locus}->{ids}->{$locusId}->{childId}->{$a} <=> $json->{locus}->{ids}->{$locusId}->{childId}->{$b} }
            keys %{ $json->{locus}->{ids}->{$locusId}->{childId} }
          )
        {
            printer_subngff( $priId, "primary", $outfh );
            foreach my $processId (
                sort { $json->{primary}->{ids}->{$priId}->{childId}->{$a} <=> $json->{primary}->{ids}->{$priId}->{childId}->{$b} }
                keys %{ $json->{primary}->{ids}->{$priId}->{childId} }
              )
            {
                if ( !exists $Done{$processId} ) {
                    printer_subngff( $processId, "processed", $outfh );
                    $Done{$processId} = 1;
                }
                foreach my $productId (
                    sort { $json->{processed}->{ids}->{$processId}->{childId}->{$a} <=> $json->{processed}->{ids}->{$processId}->{childId}->{$b} }
                    keys %{ $json->{processed}->{ids}->{$processId}->{childId} }
                  )
                {
                    if ( !exists $Done{$productId} ) {
                        printer_subngff( $productId, "product", $outfh );
                        $Done{$productId} = 1;
                    }
                }
            }
        }
    }
}

sub printer_subngff {
    my ( $id, $level, $outfh ) = @_;
    my ( $seqId, $strand, $featureType, $parentId, $start, $end, $attributes );
    #seqId	start	end	strand	level	featureType	id	parentId	attributes
    #				Chr1	6124	6578	+	primary	region	L01.pt7	L01	primary_biotype=primary_transcript;transcribed=TRUE;
    $seqId       = $json->{$level}->{"ids"}->{$id}->{"seqId"};
    $strand      = $json->{$level}->{"ids"}->{$id}->{"strand"};
    $featureType = $json->{$level}->{"ids"}->{$id}->{"featureType"};
    my @parentIds = keys %{ $json->{$level}->{"ids"}->{$id}->{parentId} };
    $parentId = join( ",", @parentIds );
    if ( exists $json->{$level}->{"ids"}->{$id}->{"regulator"} && $json->{$level}->{"ids"}->{$id}->{"regulator"} ne "" ) {
        foreach my $regulator_biotype ( keys %{ $json->{$level}->{"ids"}->{$id}->{"regulator"} } ) {
            my @Values = @{ $json->{$level}->{"ids"}->{$id}->{"regulator"}->{$regulator_biotype} };
            for my $values (@Values) {
                ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$values;
                print $outfh join( "\t", $seqId, $start, $end, $strand, $level, $featureType, $id, $parentId, $attributes ) . "\n";
            }
        }
    }
    if ( !exists $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} ) {
        print STDERR "Warning: $id has not nonRegulator value\n";
    }
    elsif ( $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} ne "" ) {
        my @Values = @{ $json->{$level}->{"ids"}->{$id}->{"nonRegulator"} };
        foreach my $values (@Values) {
            ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$values;
            if ( $level eq "locus" ) {
                $attributes = "." if !defined $attributes || $attributes eq "";
            }
            print $outfh join( "\t", $seqId, $start, $end, $strand, $level, $featureType, $id, $parentId, $attributes ) . "\n";
        }
    }
}

sub printer_gtf {
    my ( $json, $outfh ) = @_;
    my @Headers = @{ $json->{"header"} };
    for my $each (@Headers) {
        next if $each eq "";
        if ( $each =~ /genome-build|genome-version|genome-date|genome-build-accession|genebuild-last-updated/ ) {
            print $outfh "$each\n";
        }
    }
    my ( $seqId, $start, $end, $strand, $featureType, $parentId, $attributes, $gid );
    if (!defined $json->{"processed"}->{"turns"}){
	print STDERR "please add --informat parameter. exit...\n";
	exit;
    }
    my @ids = @{ $json->{"processed"}->{"turns"} };
    for my $id (@ids) {    #process id
        $seqId       = $json->{"processed"}->{"ids"}->{$id}->{"seqId"};
        $strand      = $json->{"processed"}->{"ids"}->{$id}->{"strand"};
        $featureType = $json->{"processed"}->{"ids"}->{$id}->{"featureType"};
        if ( $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator"} ne "" ) {
            $start = $json->{"processed"}->{"ids"}->{$id}->{"min_start"};
            $end   = $json->{"processed"}->{"ids"}->{$id}->{"max_end"};
            my $Attr_str;
            foreach my $tag ( sort keys %{ $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"} } ) {
                my $value = $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{$tag};
                next if $tag eq "Parent" || $tag eq "featureUid";
                next if $tag eq "transcript_id";
                $Attr_str .= " $tag \"$value\";";
            }

            if ( exists $Gene{tid2gid}{$id} ) {
                $gid = $Gene{tid2gid}{$id};
            }
            else {
                $gid = ( keys %{ $json->{processed}->{ids}->{$id}->{parentId} } )[0];

            }
            $Attr_str = "gene_id \"$gid\"; transcript_id \"$id\";" . $Attr_str;    #modified 0305

            if ( $featureType ne "polyA" ) {
                print $outfh join( "\t", $seqId, ".", "transcript", $start, $end, ".", $strand, ".", $Attr_str ) . "\n";
            }
            my @Values = @{ $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator"} };

            # exon or intron lines /lncRNA /circleRNA
            my %count;
            for my $values (@Values) {
                ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$values;
                $attributes =~ s/=/ "/g;
                $attributes =~ s/;/"; /g;
                $attributes =~ s/""/"/g;
                $count{$featureType}++;
                my $featureNum = $count{$featureType};
                if ( $attributes !~ /\s${featureType}_number\s/ ) {
                    $attributes .= "${featureType}_number \"$featureNum\";";
                }
                if ( $attributes !~ /\stranscript_id\s/ ) {
                    $attributes = $Attr_str . $attributes;    #modified 0305
                }
                else {                                        #modified 0305
                    $attributes =~ s/ transcript_id \".*?\"//;
                    $attributes = "gene_id \"$gid\"; transcript_id \"$id\"; " . $attributes;    #modified 0305
                }    #modified 0305
                print $outfh join( "\t", $seqId, ".", $featureType, $start, $end, ".", $strand, ".", $attributes ) . "\n";
            }

        }

    }
}

my $hashSort = {};

sub sorter {
    foreach my $level ( keys %{$json} ) {
        next if $level eq "header" || $level eq "input_format";
        if ( exists $json->{$level}->{ids} ) {
            $hashSort = $json->{$level}->{ids};
            foreach my $id ( sort keys %{$hashSort} ) {
                if ( !exists $hashSort->{$id}->{seqId} ) {
                    delete $hashSort->{$id};
                }
            }
            my @Turns = sort by_seq_start_end keys %{$hashSort};
            @{ $json->{$level}->{"turns"} } = @Turns;
        }
    }
}

sub by_seq_start_end {
    if    ( $hashSort->{$a}->{"seqId"} lt $hashSort->{$b}->{"seqId"} ) { -1 }
    elsif ( $hashSort->{$a}->{"seqId"} gt $hashSort->{$b}->{"seqId"} ) { 1 }
    else {
        if    ( $hashSort->{$a}->{"strand"} lt $hashSort->{$b}->{"strand"} ) { -1 }
        elsif ( $hashSort->{$a}->{"strand"} gt $hashSort->{$b}->{"strand"} ) { 1 }
        else {
            if    ( $hashSort->{$a}->{"min_start"} < $hashSort->{$b}->{"min_start"} ) { -1 }
            elsif ( $hashSort->{$a}->{"min_start"} > $hashSort->{$b}->{"min_start"} ) { 1 }
            else {
                if    ( $hashSort->{$a}->{"max_end"} < $hashSort->{$b}->{"max_end"} || $hashSort->{$a}->{"max_end"} == $hashSort->{$b}->{"max_end"} ) { -1 }
                elsif ( $hashSort->{$a}->{"max_end"} > $hashSort->{$b}->{"max_end"} )                                                                 { 1 }
            }
        }
    }
}

sub printer_gff3 {
    my ( $json, $outfh ) = @_;
    CDSphase();
    my @Headers = @{ $json->{"header"} };
    for my $each (@Headers) {
        next if $each eq "";
        if ( $each =~ /genome-build|genome-version|genome-date|genome-build-accession|genebuild-last-updated/ ) {
            print $outfh "$each\n";
        }
    }
    if ( exists $json->{"processed"} ) {
        my @ids = @{ $json->{"processed"}->{"turns"} };
        subprinter_gff3( "novalue", $json, $outfh, @ids );
    }
}

sub subprinter_gff3 {
    my $gene  = shift @_;
    my $json  = shift @_;
    my $outfh = shift @_;
    my @ids   = @_;
    my ( $seqId, $start, $end, $strand, $featureType, $featureType_out, $parentId, $attributes );
    for my $id (@ids) {    #process id
        $seqId       = $json->{"processed"}->{"ids"}->{$id}->{"seqId"};
        $strand      = $json->{"processed"}->{"ids"}->{$id}->{"strand"};
        $featureType = $json->{"processed"}->{"ids"}->{$id}->{"featureType"};
        if ( $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator"} ne "" ) {
            $start = $json->{"processed"}->{"ids"}->{$id}->{"min_start"};
            $end   = $json->{"processed"}->{"ids"}->{$id}->{"max_end"};
            my $Attr_str;
            my $Attr_str_gene = "";
            foreach my $tag ( sort keys %{ $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"} } ) {
                my $value = $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{$tag};
                next if $tag eq "featureUid";
                next if $tag eq "Parent";
                $Attr_str .= "$tag=$value;";
                if ( $tag =~ /gene/i ) {
                    $Attr_str_gene .= "$tag=$value;";
                }
            }

            #gene line if gtf2gff3
            my ( $gid, $start, $end );
            if ( $gene ne "novalue" ) {
                $gid   = $gene;
                $start = $Gene{gid2productid}{$gene}{min_start};
                $end   = $Gene{gid2productid}{$gene}{max_end};
            }
            else {
                if ( !exists $Gene{tid2gid}{$id} ) {
                    $gid = ( keys %{ $json->{processed}->{ids}->{$id}->{parentId} } )[0];
                }
                else {
                    $gid = $Gene{tid2gid}{$id};
                }
                $start = $Gene{gidstart}{$gid};
                $end   = $Gene{gidend}{$gid};
            }
            if ( !exists $Gene{gidDone}{$gid} ) {
                print $outfh join( "\t", $seqId, ".", "gene", $start, $end, ".", $strand, ".", "ID=$gid;$Attr_str_gene" ) . "\n";
            }
            $Gene{gidDone}{$gid} = 1;

            if ( exists $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype}
                && $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype} eq "miRNA" )
            {
                $featureType_out = "primary_transcript";
            }
            elsif ( exists $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype}
                && $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype} eq "lncRNA" )
            {
                $featureType_out = "lncRNA";
            }
            elsif ( exists $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype}
                && $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype} eq "circleRNA" )
            {
                $featureType_out = "circleRNA";
            }
            else {
                $featureType_out = "mRNA";
            }
            $start = $json->{"processed"}->{"ids"}->{$id}->{"min_start"};    #modified 0305  ngff ->gff3 mRNA is wrong
            $end   = $json->{"processed"}->{"ids"}->{$id}->{"max_end"};      #modified 0305  ngff ->gff3
            if ( $gene ne "novalue" ) {
                print $outfh join( "\t", $seqId, ".", $featureType_out, $start, $end, ".", $strand, ".", "ID=$id;Parent=$gene;$Attr_str" ) . "\n";
            }
            else {
                my $gid;
                if ( !exists $Gene{tid2gid}{$id} ) {
                    $gid = ( keys %{ $json->{processed}->{ids}->{$id}->{parentId} } )[0];
                }
                else {
                    $gid = $Gene{tid2gid}{$id};
                }
                print $outfh join( "\t", $seqId, ".", $featureType_out, $start, $end, ".", $strand, ".", "ID=$id;Parent=$gid;$Attr_str" ) . "\n";
            }

            my @Values = @{ $json->{"processed"}->{"ids"}->{$id}->{"nonRegulator"} };

            # exon or intron lines /lncRNA /circleRNA
            my %count;
            for my $values (@Values) {
                ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$values;
                $count{$featureType}++;
                my $featureNum = $count{$featureType};
                if ( $attributes =~ /${featureType}_number/ ) {
                    print $outfh
                      join( "\t", $seqId, ".", $featureType, $start, $end, ".", $strand, ".", "ID=${id}${featureType}${featureNum};Parent=$id;${attributes}" )
                      . "\n";
                }
                else {
                    print $outfh join( "\t",
                        $seqId, ".", $featureType, $start, $end, ".", $strand, ".",
                        "ID=${id}${featureType}${featureNum};Parent=$id;${attributes}${featureType}_number=$featureNum;" )
                      . "\n";
                }
            }
            foreach my $productid ( sort keys %{ $json->{"processed"}->{"ids"}->{$id}->{"childId"} } ) {
                my @Values = @{ $json->{"product"}->{"ids"}->{$productid}->{"nonRegulator"} };
                for my $num ( 0 .. $#Values ) {
                    my $values = $Values[$num];
                    ( $seqId, $start, $end, $strand, $featureType, $attributes ) = @$values;
                    next if ( $start < 0 && $end < 0 );
                    next if ( $featureType eq "regulator" );
                    my $phase = ".";
                    if ( $featureType eq "CDS" ) {
                        $phase = $Gene{CDSphase}{$productid}{$num};
                    }

                    elsif ( exists $json->{"product"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype}
                        && $json->{"product"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{gene_biotype} eq "miRNA" )
                    {
                        $featureType = "miRNA";
                    }
                    if ( exists $json->{"product"}->{"ids"}->{$id}->{"nonRegulator_shared_attributes"}->{protein_id} ) {
                        $attributes = "protein_id=$productid;" . $attributes;
                    }
                    print $outfh join( "\t", $seqId, ".", $featureType, $start, $end, ".", $strand, $phase, "ID=${productid};Parent=$id;$attributes" ) . "\n";
                }

            }

        }

    }

}

sub printer_json { }
sub printer_ntsv { }

sub convert_AbsPos2RelPos {
    my ( $pos, $abs_start, $rel_start, $rel_end ) = @_;
    my $rel_pos = $rel_start + ( $pos - $abs_start ) * ( $rel_end - $rel_start ) / abs( $rel_end - $rel_start );
    return $rel_pos;
}

sub convert_RelPos2AbsPos {
    my ( $pos, $abs_start, $rel_start ) = @_;
    my $abs_pos = $rel_start - $pos + $abs_start;
    return $abs_pos;
}

sub Uniq_array {
    my @array = @_;
    my %count;
    my @uniq_order = grep { ++$count{$_} < 2; } @array;
    return @uniq_order;
}

sub ArrayToHash {
    my $str = shift @_;
    my %hash;
    my @tmp = split /;|=/, $str;
    for ( my $i = 0 ; $i < $#tmp ; $i += 2 ) {
        my $tag   = $tmp[$i];
        my $value = $tmp[ $i + 1 ];
        $hash{$tag} = $value;
    }
    return %hash;
}

1;
