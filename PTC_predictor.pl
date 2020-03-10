#requirements
use strict; 
use warnings;
no warnings ('uninitialized', 'substr');
use Bio::EnsEMBL::Registry;

# get input from terminal
my ($rownumber, $field) = @ARGV;

#load ensembl
my $registry = "Bio::EnsEMBL::Registry";
$registry->load_registry_from_db(
    -host => 'ensembldb.ensembl.org',
    -user => 'anonymous',
    -species => 'homo_sapiens',
    );
my $slice_adaptor = $registry->get_adaptor('Human', 'Core', 'Slice');

#load input                   
my ($locus, $muttype, $mutsubtype, $mutchange) = split /-/, $field;
my ($chr, $pos) = split /:/, $locus;
my $chr_num = substr $chr, 3, (length($chr) - 3);

# get info from ensembl
my $slice = $slice_adaptor->fetch_by_region('chromosome', $chr_num, $pos, $pos);
my @genes = $slice->get_all_Genes();
my @transcripts_locus = $slice->get_all_Transcripts();

my @transcripts_locus_ids;
foreach my $transcript_locus ( @{ $slice->get_all_Transcripts() }){
    push @transcripts_locus_ids, $transcript_locus->stable_id();
}

foreach my $gene ( @{ $slice->get_all_Genes()} ) {
    
    foreach my $transcript ( @{ $gene->get_all_Transcripts()} ) {
        my $transcript_id = $transcript->stable_id();

        # transcripts for each gene that correspond to locus
        if ( grep { $_ eq $transcript_id} @transcripts_locus_ids) {
            my $transcript_biotype = $transcript->biotype();  

            # set per transcript output to zero
            my $out_trans_chr = $chr_num;
            my $out_trans_pos = $pos;
            my $out_trans_mutsubtype = $mutsubtype;
            my $out_trans_muttype = $muttype;
            my $out_trans_mutchange = $mutchange;
            my $out_trans_gene = $transcript->get_Gene()->external_name();
            my $out_trans_gene_id = $transcript->get_Gene()->display_id();
            my $out_trans_transcript_id = $transcript_id;
            my $out_trans_biotype = $transcript_biotype;
            my $out_trans_strand = $transcript->strand();
            my $out_trans_n_exons = "undef";
            my $out_trans_n_coding_exons = "undef";
            my $out_trans_mut_exon = "undef";
            my $out_trans_mut_n_exon = "undef";
            my $out_trans_mut_n_coding_exon = "undef";
            my $out_trans_mut_exon_start = "undef";
            my $out_trans_mut_exon_end = "undef";
            my $out_trans_mut_exon_seq = "undef";
            my $out_trans_microsat_length = "undef";
            my $out_trans_microsat_allel = "undef";
            my $out_trans_ptc_found = "undef";
            my $out_trans_ptc_seq = "undef";
            my $out_trans_ptc_distance = "undef";
            my $out_trans_ptc_exon = "undef";
            my $out_trans_ptc_n_exon = "undef";
            my $out_trans_ptc_n_coding_exon = "undef";
            my $out_trans_ptc_exon_start = "undef";
            my $out_trans_ptc_exon_end = "undef";
            my $out_trans_ptc_exon_seq = "undef";                               

            # PTC analysis on selected transcripts                         
            if ($mutsubtype eq "Frameshift" && $muttype ne "NONE" && $muttype ne "UNDEF"){
                # PTC analysis on frameshift transcripts

                # check exon containing mut
                my @exons =  @{ $transcript->get_all_Exons() };
                my $n_exons = scalar(@exons);
                $out_trans_n_exons = $n_exons;
                my $exon_count = 0;
                my $coding_exon_count = 0;
                
                foreach my $exon (@exons){
                    
                    # check if mutation is in exon
                    $exon_count++;
                        
                    if ($exon->is_coding($transcript) == 1){
                        $coding_exon_count++;
                        my $exon_start = $exon->seq_region_start();
                        my $exon_end = $exon->seq_region_end();

                        # check if mutation is in coding part exon
                        if ($exon_start <= $pos && $exon_end >= $pos){
                            my $exon_strand = $exon->seq_region_strand();
                            $out_trans_mut_exon = $exon->stable_id();
                            $out_trans_mut_n_exon = $exon_count;
                            $out_trans_mut_n_coding_exon = $coding_exon_count;
                            my $distance_from_coding_start = $pos - $exon_start;
                            my $distance_from_coding_end = $exon_end - $pos;
                            $out_trans_mut_exon_start = $distance_from_coding_start;
                            $out_trans_mut_exon_end = $distance_from_coding_end;
                            my $exon_seq = $exon->seq()->seq();
                            $out_trans_mut_exon_seq  = $exon_seq;
                            
                            if($exon_strand == 1){
                                
                                if($coding_exon_count == 1){
                    
                                    my $coding_part = $exon->coding_region_start($transcript) + $distance_from_coding_start;
                                    if ($coding_part >= 0){
                                        
                                    } 
                                }
                                
                                my $mut_seq = substr $exon_seq, ($distance_from_coding_start + 1), $mutchange;
                                my $mut_seq_flanked_pos = substr $exon_seq, ($distance_from_coding_start + 1), 10;
                                
                                $out_trans_microsat_allel = $mut_seq;

                                # calculate microsattelite repeats
                                my $micorsat_repeats = 1;
                                my $checkpoint = $distance_from_coding_start + 1 + $mutchange;
                                my $microsat_check_seq = substr $exon_seq, $checkpoint, $mutchange;
                                while($microsat_check_seq eq $mut_seq && length($microsat_check_seq) > 0){
                                    $micorsat_repeats++;
                                    $checkpoint = $checkpoint + $mutchange;
                                    $microsat_check_seq = substr $exon_seq, $checkpoint, $mutchange;
                                }
                                
                                $out_trans_microsat_length = $micorsat_repeats;

                                # Calculate distance to PTC
                                my $PTC_search_checkpoint;
                                my $exon_phase = $exon->phase();
                                my $nucl_offset;

                                if ($exon_phase != -1){
                                    $nucl_offset = 3 - $exon_phase;
                                }else{
                                    $nucl_offset = (($exon->coding_region_start($transcript) * -1) -1) % 3;
                                }
                                my $frame_distance = $distance_from_coding_start + 1 - $nucl_offset;
                                my $remainder = $frame_distance % 3;
                                
                                if($muttype eq "DEL"){
                                    $PTC_search_checkpoint = $distance_from_coding_start + 1 + $mutchange + (3 - $remainder);
                                }elsif($muttype eq "INS"){
                                    $PTC_search_checkpoint = $distance_from_coding_start + 1 - $mutchange + (3 - $remainder);
                                }else{
                                }
                                my $PTC_search = substr $exon_seq, $PTC_search_checkpoint, 3;
                                my $AA_count = 0;
                                until($PTC_search eq "TAG" || $PTC_search eq "TAA" || $PTC_search eq "TGA" || length($PTC_search) < 3){
                                    $AA_count++;
                                    $PTC_search_checkpoint = $PTC_search_checkpoint + 3;
                                    $PTC_search = substr $exon_seq, $PTC_search_checkpoint, 3;
                                    
                                }
                                # check if PTC is in current exon or continue to next exon
                                if($PTC_search eq "TAG" || $PTC_search eq "TAA" || $PTC_search eq "TGA"){
                                    $out_trans_ptc_found = "TRUE";
                                    $out_trans_ptc_seq = $PTC_search;
                                    $out_trans_ptc_distance = $AA_count;
                                    $out_trans_ptc_exon = $exon->stable_id();
                                    $out_trans_ptc_n_exon = $exon_count;
                                    $out_trans_ptc_n_coding_exon = $coding_exon_count;
                                    $out_trans_ptc_exon_start = $distance_from_coding_start;
                                    $out_trans_ptc_exon_end = $distance_from_coding_end;
                                    $out_trans_ptc_exon_seq = $exon->seq()->seq();
                                }else{
                                    my $exon_to_check = $coding_exon_count;
                                    my $exons_checked = 0;
                                    my $AA_distance = $AA_count;
                                    my $stop_codon_found = "FALSE";
                                    my $last_coding_exon_checked = "FALSE";
                                    my $next_PTC_search_offset = length($PTC_search);

                                    until($stop_codon_found eq "TRUE" || $last_coding_exon_checked eq "TRUE"){
                                        $exon_to_check++;
                                        $exons_checked++;
                                        my $next_coding_exon_counter = 0;
                                        my $next_AA_count = 0;
                                        foreach my $next_exon (@exons){
                                            
                                            if ($next_exon->is_coding($transcript) == 1){
                                                $next_coding_exon_counter++;
                                                
                                                if($next_coding_exon_counter == $exon_to_check){
                                                
                                                    my $next_exon_seq = $next_exon->seq()->seq();
                                                    my $next_PTC_checkpoint = 3 - $next_PTC_search_offset;
                                                    my $next_PTC_search = substr $next_exon_seq, $next_PTC_checkpoint, 3;
                                                    
                                                    until($next_PTC_search eq "TAG" || $next_PTC_search eq "TAA" || $next_PTC_search eq "TGA" || length($next_PTC_search) < 3){
                                                        $next_AA_count++;
                                                        $AA_distance++;
                                                        $next_PTC_checkpoint = $next_PTC_checkpoint + 3;
                                                        $next_PTC_search = substr $next_exon_seq, $next_PTC_checkpoint, 3;
                                                    }
                                                    if($next_PTC_search eq "TAG" || $next_PTC_search eq "TAA" || $next_PTC_search eq "TGA"){
                                                        $stop_codon_found = "TRUE";
                                                        $out_trans_ptc_found = "TRUE";
                                                        $out_trans_ptc_seq = $next_PTC_search;
                                                        $out_trans_ptc_distance = $AA_distance;
                                                        $out_trans_ptc_exon = $next_exon->stable_id();
                                                        $out_trans_ptc_n_exon = $exon_count + $exons_checked;
                                                        $out_trans_ptc_n_coding_exon = $coding_exon_count + $exons_checked;
                                                        $out_trans_ptc_exon_start = $next_PTC_checkpoint;
                                                        $out_trans_ptc_exon_end = $next_exon->length() - $next_PTC_checkpoint;
                                                        $out_trans_ptc_exon_seq = $next_exon->seq()->seq();
                                                    }
                                                    $next_PTC_search_offset = length($next_PTC_search);
                                                }

                                            }
                                        }
                                        
                                        if($next_coding_exon_counter < $exon_to_check){
                                            $last_coding_exon_checked = "TRUE";
                                            $out_trans_ptc_found = "FALSE";
                                            $out_trans_ptc_distance = $AA_distance;
                                        }
                                    }

                                }

                            }elsif($exon_strand == -1){
                                my $mut_seq = substr $exon_seq, ($distance_from_coding_end - 1), $mutchange;
                                my $mut_seq_flanked_pos = substr $exon_seq, ($distance_from_coding_end - 1 ), 10;
                                $out_trans_microsat_allel = $mut_seq;

                                # calculate microsattelite repeats
                                my $micorsat_repeats = 1;
                                my $checkpoint = $distance_from_coding_end - 1 - $mutchange;
                                my $microsat_check_seq = substr $exon_seq, $checkpoint, $mutchange;
                                while($microsat_check_seq eq $mut_seq && length($microsat_check_seq) > 0){
                                    $micorsat_repeats++;
                                    $checkpoint = $checkpoint - $mutchange;
                                    $microsat_check_seq = substr $exon_seq, $checkpoint, $mutchange;
                                }
                                # print $micorsat_repeats, "\t";
                                $out_trans_microsat_length = $micorsat_repeats;

                                # Calculate distance to PTC
                                my $PTC_search_checkpoint;
                                my $exon_phase = $exon->phase();
                                my $nucl_offset;

                                if ($exon_phase != -1){
                                    $nucl_offset = 3 - $exon_phase;
                                }else{
                                    $nucl_offset = (($exon->coding_region_start($transcript) * -1) -1) % 3;
                                }
                                my $frame_distance = $distance_from_coding_end - 1 - $nucl_offset;
                                my $remainder = $frame_distance % 3;
                                
                                if($muttype eq "DEL"){
                                    $PTC_search_checkpoint = $distance_from_coding_end - 1 + $mutchange + (3 - $remainder);
                                }elsif($muttype eq "INS"){
                                    $PTC_search_checkpoint = $distance_from_coding_end - 1 - $mutchange + (3 - $remainder);
                                }else{
                                }
                                my $PTC_search = substr $exon_seq, $PTC_search_checkpoint, 3;
                                my $AA_count = 0;
                                until($PTC_search eq "TAG" || $PTC_search eq "TAA" || $PTC_search eq "TGA" || length($PTC_search) < 3){
                                    $AA_count++;
                                    $PTC_search_checkpoint = $PTC_search_checkpoint + 3;
                                    $PTC_search = substr $exon_seq, $PTC_search_checkpoint, 3;
                                }
                                
                                if($PTC_search eq "TAG" || $PTC_search eq "TAA" || $PTC_search eq "TGA"){
                                    $out_trans_ptc_found = "TRUE";
                                    $out_trans_ptc_seq = $PTC_search;
                                    $out_trans_ptc_distance = $AA_count;
                                    $out_trans_ptc_exon = $exon->stable_id();
                                    $out_trans_ptc_n_exon = $exon_count;
                                    $out_trans_ptc_n_coding_exon = $coding_exon_count;
                                    $out_trans_ptc_exon_start = $distance_from_coding_start;
                                    $out_trans_ptc_exon_end = $distance_from_coding_end;
                                    $out_trans_ptc_exon_seq = $exon->seq()->seq();
                                }else{
                                    my $exon_to_check = $coding_exon_count;
                                    my $exons_checked = 0;
                                    my $stop_codon_found = "FALSE";
                                    my $last_coding_exon_checked = "FALSE";
                                    my $next_PTC_search_offset = length($PTC_search);

                                    until($stop_codon_found eq "TRUE" || $last_coding_exon_checked eq "TRUE"){
                                        $exon_to_check++;
                                        $exons_checked++;
                                        my $AA_distance = $AA_count;
                                        my $next_coding_exon_counter = 0;
                                        my $next_AA_count = 0;
                                        foreach my $next_exon (@exons){
                                            
                                            if ($next_exon->is_coding($transcript) == 1){
                                                $next_coding_exon_counter++;
                                                
                                                if($next_coding_exon_counter == $exon_to_check){
                                                    my $next_exon_seq = $next_exon->seq()->seq();
                                                    my $next_PTC_checkpoint = 3 - $next_PTC_search_offset;
                                                    my $next_PTC_search = substr $next_exon_seq, $next_PTC_checkpoint, 3;
                                                    
                                                    until($next_PTC_search eq "TAG" || $next_PTC_search eq "TAA" || $next_PTC_search eq "TGA" || length($next_PTC_search) < 3){
                                                        $next_AA_count++;
                                                        $AA_distance++;
                                                        $next_PTC_checkpoint = $next_PTC_checkpoint + 3;
                                                        $next_PTC_search = substr $next_exon_seq, $next_PTC_checkpoint, 3;
                                                    }
                                                    if($next_PTC_search eq "TAG" || $next_PTC_search eq "TAA" || $next_PTC_search eq "TGA"){
                                                        $stop_codon_found = "TRUE";
                                                        $out_trans_ptc_found = "TRUE";
                                                        $out_trans_ptc_seq = $next_PTC_search;
                                                        $out_trans_ptc_distance = $AA_distance;
                                                        $out_trans_ptc_exon = $next_exon->stable_id();
                                                        $out_trans_ptc_n_exon = $exon_count + $exons_checked;
                                                        $out_trans_ptc_n_coding_exon = $coding_exon_count + $exons_checked;
                                                        $out_trans_ptc_exon_start = $next_PTC_checkpoint;
                                                        $out_trans_ptc_exon_end = $next_exon->length() - $next_PTC_checkpoint;
                                                        $out_trans_ptc_exon_seq = $next_exon->seq()->seq();
                                                    }
                                                    $next_PTC_search_offset = length($next_PTC_search);
                                                }

                                            }
                                        }
                                        
                                        if($next_coding_exon_counter < $exon_to_check){
                                            $last_coding_exon_checked = "TRUE";
                                            $out_trans_ptc_found = "FALSE";
                                            $out_trans_ptc_distance = $AA_distance;
                                        }
                                    }

                                }
                            }else{

                            }

                        }

                    }
                    $out_trans_n_coding_exons = $coding_exon_count
                }
                
            }elsif ($muttype eq "SNP"){
                ### microsat analysis on SNPs
                my $micorsat_repeats = 0;
                my $check_pos = $pos + 1;
                my $initial_slice = $slice_adaptor->fetch_by_region('chromosome', $chr_num, $check_pos, $check_pos);
                my $initial_slice_seq = $initial_slice->seq();
                my $check_slice = $initial_slice;
                my $check_slice_seq = $check_slice->seq();
                while($initial_slice_seq eq $check_slice_seq){
                    $micorsat_repeats++;
                    $check_pos++;
                    $check_slice = $slice_adaptor->fetch_by_region('chromosome', $chr_num, $check_pos, $check_pos);
                    $check_slice_seq = $check_slice->seq();                           
                }
                $out_trans_microsat_length = $micorsat_repeats;
                $out_trans_microsat_allel = $check_slice_seq;
            }else{
               
            }

            #print output per transcript.
            print $out_trans_chr, "\t";
            print $out_trans_pos, "\t";
            print $out_trans_mutsubtype, "\t";
            print $out_trans_muttype, "\t";
            print $out_trans_mutchange, "\t";
            print $out_trans_gene, "\t";
            print $out_trans_gene_id, "\t";
            print $out_trans_transcript_id, "\t";
            print $out_trans_biotype, "\t";
            print $out_trans_strand, "\t";
            print $out_trans_n_exons, "\t";
            print $out_trans_n_coding_exons, "\t";
            print $out_trans_mut_exon, "\t";
            print $out_trans_mut_n_exon, "\t";
            print $out_trans_mut_n_coding_exon, "\t";
            print $out_trans_mut_exon_start, "\t";
            print $out_trans_mut_exon_end, "\t";
            #print $out_trans_mut_exon_seq, "\t";
            print $out_trans_microsat_length, "\t";
            print $out_trans_microsat_allel, "\t";
            print $out_trans_ptc_found, "\t";
            print $out_trans_ptc_seq, "\t";
            print $out_trans_ptc_distance, "\t";
            print $out_trans_ptc_exon, "\t";
            print $out_trans_ptc_n_exon, "\t";
            print $out_trans_ptc_n_coding_exon, "\t";
            print $out_trans_ptc_exon_start, "\t";
            print $out_trans_ptc_exon_end, "\n";
            # print $out_trans_ptc_exon_seq, "\n";
        }
        
    }
    
} 