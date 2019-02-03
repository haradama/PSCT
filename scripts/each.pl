use strict;
use G;
use SUBHAROUTINE;

my $dir = ".";
my $dir = $ARGV[0];
my $data = $ARGV[1]; # 'A0' 'R0' 'R4'
my $del_key = '[^CDEFHKNQYacgt]'; # 2-fold degenerate amino acids
if($data =~ /R0|R4/){ $del_key = '[^ACDEFGHIKLNPQRSTVYacgt]'; }
elsif($data =~ /A0/){ $del_key = '[^ACDEFGHIKLMNPQRSTVWYacgtU]'; } # U Sec Selenocysteine
my $off = '0nltk'; # '200nltk'

##### analyses
opendir(DIR, $dir) || die "directory open error:$!";
foreach my $file (sort readdir(DIR)){
    next if($file !~ /(.+)\.(gbk|gbff)$/); my $filename = $1;
    #next if($file !~ /NC_010935/); # test # plasmid pCNB
    #next if($file !~ /SRR1749435/); # test # MetaSUB
    my $gb = load("$dir/$file","no msg");
    $filename =~ s/_prokka_//;
    $gb->{haruo}->{filename} = $filename;
    $gb->{ACCESSION} = $filename if($gb->{ACCESSION} =~ /^VERSION$/);
    print "\n$filename: $gb->{DEFINITION}\n";
    
    &subharoutine::makedata($gb, "$dir/$filename.gbk", 0); # $oriter = 1
    my $parameter = ['Laa','aroma','gravy','mmw','gcc3','gtc3','Hgc3','Ew','P2','cai']; #  'enc','icdi',,'gcs3','CMHFYW', 'fop','Lc',  'BgC','BgH','E_g','phx','pa','A_g']; 'a_c3','c_c3','g_c3','t_c3', 'gcc',
    $parameter = [@$parameter, 'start','end','direction', 'locus_tag','protein_id','gene','product', 'highlyExpressed', 'leading'];  # 'db_xref','GI','GeneID','EC_number', 'cog04','cog25','COG','COGproduct', 'High', 'tag','ACCESSION', 
    
    &codon_compiler($gb, -output=>'stdout', -del_key=>$del_key, -data=>$data, -startcodon=>0, -stopcodon=>0);

    print "\n ndc = ", $gb->{All}->{ndc}, "\n"; next if($gb->{All}->{ndc} > 59);

    foreach my $cds ($gb->cds()){
	$gb->{$cds}->{ACCESSION} = $gb->{ACCESSION};
	if($gb->{$cds}->{Laa} < $off){ $gb->{$cds}->{on} = 0; next; }
	&codon_compiler($gb, -output=>'n', -del_key=>$del_key, -data=>$data, -startcodon=>0, -stopcodon=>0, -id=>$cds);
	if($off =~ /nltk/ && $gb->{$cds}->{nltk}){ $gb->{$cds}->{on} = 0; next; }
	&aaui($gb, -output=>'n', -id=>$cds);
	#&bui($gb, -output=>'n', -position=>'', -del_key=>$del_key, -id=>$cds);
	&bui($gb, -output=>'n', -position=>3, -del_key=>$del_key, -id=>$cds);
    }
    &Ew($gb, -output=>'n');
    &P2($gb, -output=>'n');
    #my $w_val = &w_value($gb, -output=>'n', -sharp=>1); 
    my $w_val = &w_value($gb, -output=>'n', -sharp=>1, -tag=>"KOGproduct"); # for paramecium
    &cai($gb, -output=>'n', -w_values=>$w_val);
    #&cai($gb, -output=>'n', -w_output=>'n');
    #&cai($gb, -output=>'n', -w_output=>'n', -tai=>1);
    #&phx($gb, -output=>'n'); # E_g
    
    #open(OUT, ">$filename.$data.txt");
    open(OUT, ">$filename\_$off$data.txt");
    my @keyAll = sort keys %{$gb->{All}->{$data}};
    print OUT join("\t", (@keyAll, @$parameter) ),"\n";
    foreach my $cds ($gb->cds()){
	if($cds =~ /FEATURE([0-9]+)/){ $gb->{$cds}->{tag} = $filename."_".sprintf("%04d",$1); }
	#$gb->{$cds}->{gene} =~ s/,/_/g; $gb->{$cds}->{product} =~ s/,/_/g; $gb->{$cds}->{function} =~ s/,/_/g; $gb->{$cds}->{note} =~ s/,/_/g; $gb->{$cds}->{COGproduct} =~ s/,/_/g;
	my @tmp;
	foreach (@keyAll){ push(@tmp, $gb->{$cds}->{$data}->{$_} || 0); }
	foreach (@$parameter){ push(@tmp, $gb->{$cds}->{$_}); }
	print OUT join("\t", @tmp), "\n";
    }
	close(OUT);
}
closedir(DIR);

__END__
http://kobesearch.cpan.org/htdocs/Bio-Glite/Bio/Glite.html

data="R0"
#data="A0"
dir="."; 
(time perl -I/$HOME/g-language-1.9.1/lib ~/scripts/each.pl $dir $data &) >& log.each.$data.txt
