my @site = &seq2mismatch($gRNA_seq,$ref_seq);
my $gRNA_score = &sites2score(@site );



sub seq2mismatch{
	my ($query_dna,$ref_dna) = @_;
	my @out = ();
	# end with letter
	# end with number
	
	if(length($query_dna) != length($ref_dna)){
		die "not same length\n";
	}
	
	for (my $i=0;$i<length($query_dna);$i++){
		if((substr $query_dna,$i,1) ne (substr $ref_dna,$i,1) ){
			if($i<20){
				push @out,$i;
			}			
		}
	}
	
	return @out;
}


sub sites2score{
    	my @mis_site = @_; # base 0
      my @weight_all = (0,0,0.014,0,0,0.395,0.317,0,0.389,0.079,0.445,0.508,0.613,0.851,0.732,0.828,0.615,0.804,0.685,0.583);


	          my $mis_site_num = @mis_site + 0;
	          ## score 3
	          my $score_p3 = 1;
	          if($mis_site_num>=1){
	          	  $score_p3 = 1/($mis_site_num*$mis_site_num);
	          }
	          ## score 2
	           my $score_p2 = 1;	          
	           if($mis_site_num>=2){
	           	   # others equal 1
	           	   my @distance = ();
	          	   for(my $i=0;$i<($mis_site_num-1);$i++){
	          	   	   my $dis = $mis_site[$i+1]-$mis_site[$i];
	          	   	   push @distance,$dis;
	          	   }
	          	   my $sum_dis = 0;
				   for( @distance){
					    $sum_dis = $sum_dis + $_;
				   }
				   my $mean_dis = $sum_dis/(@distance+0);
			       $score_p2 = 1/((((19-$mean_dis)/19)*4)+1);
	           }
	           
	           ## score 1
	           my $score_p1 = 1;
			   ## get score part 1
			   for(@mis_site){		
						$score_p1 = $score_p1*(1-$weight_all[$_]);
			   }
	           	           
	           my $score_hit = $score_p1*$score_p2 *$score_p3*100;
	           return $score_hit;
}# get score
    
    
    