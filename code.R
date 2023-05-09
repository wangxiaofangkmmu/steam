####################################################################################################################################
## down TCGA data case：TCGA LUAD		file：RNA-seq FPKM       Add all files to cart      metadata,download-cart
#!/usr/bin/perl -w
use strict;
use warnings;

my $file=$ARGV[0];

#use Data::Dumper;
use JSON;

my $json = new JSON;
my $js;

my %hash=();
my @normalSamples=();
my @tumorSamples=();

open JFILE, "$file";
while(<JFILE>) {
  $js .= "$_";
}
my $obj = $json->decode($js);
for my $i(@{$obj})
{
  my $file_name=$i->{'file_name'};
  my $file_id=$i->{'file_id'};
  my $entity_submitter_id=$i->{'associated_entities'}->[0]->{'entity_submitter_id'};
  $file_name=~s/\.gz//g;
  if(-f $file_name)
  {
    my @idArr=split(/\-/,$entity_submitter_id);
    if($idArr[3]=~/^0/)
    {
      push(@tumorSamples,$entity_submitter_id);
    }
    else
    {
      push(@normalSamples,$entity_submitter_id);
    }        	
    open(RF,"$file_name") or die $!;
    while(my $line=<RF>)
    {
      next if($line=~/^\n/);
      next if($line=~/^\_/);
      chomp($line);
      my @arr=split(/\t/,$line);
      ${$hash{$arr[0]}}{$entity_submitter_id}=$arr[1];
    }
    close(RF);
  }
}
#print Dumper $obj

open(WF,">mRNAmatrix.txt") or die $!;
my $normalCount=$#normalSamples+1;
  my $tumorCount=$#tumorSamples+1;
  print "normal count: $normalCount\n";
print "tumor count: $tumorCount\n";
print WF "id\t" . join("\t",@normalSamples);
print WF "\t" . join("\t",@tumorSamples) . "\n";
foreach my $key(keys %hash)
{
  print WF $key;
  foreach my $normal(@normalSamples)
  {
    print WF "\t" . ${$hash{$key}}{$normal};
  }
  foreach my $tumor(@tumorSamples)
  {
    print WF "\t" . ${$hash{$key}}{$tumor};
  }
  print WF "\n";
}
close(WF);
## 
#!/usr/bin/perl -w
use strict;
use warnings;

my $file=$ARGV[0];

#use Data::Dumper;
use JSON;

my $json = new JSON;
my $js;

my %hash=();
my @normalSamples=();
my @tumorSamples=();

open JFILE, "$file";
while(<JFILE>) {
  $js .= "$_";
}
my $obj = $json->decode($js);
for my $i(@{$obj})
{
  my $file_name=$i->{'file_name'};
  my $file_id=$i->{'file_id'};
  my $entity_submitter_id=$i->{'associated_entities'}->[0]->{'entity_submitter_id'};
  $file_name=~s/\.gz//g;
  if(-f $file_name)
  {
    my @idArr=split(/\-/,$entity_submitter_id);
    if($idArr[3]=~/^0/)
    {
      push(@tumorSamples,$entity_submitter_id);
    }
    else
    {
      push(@normalSamples,$entity_submitter_id);
    }        	
    open(RF,"$file_name") or die $!;
    while(my $line=<RF>)
    {
      next if($line=~/^\n/);
      next if($line=~/^\_/);
      chomp($line);
      my @arr=split(/\t/,$line);
      ${$hash{$arr[0]}}{$entity_submitter_id}=$arr[1];
    }
    close(RF);
  }
}
#print Dumper $obj
open(WF,">mRNAmatrix.txt") or die $!;
my $normalCount=$#normalSamples+1;
  my $tumorCount=$#tumorSamples+1;
  print "normal count: $normalCount\n";
print "tumor count: $tumorCount\n";
print WF "id\t" . join("\t",@normalSamples);
print WF "\t" . join("\t",@tumorSamples) . "\n";
foreach my $key(keys %hash)
{
  print WF $key;
  foreach my $normal(@normalSamples)
  {
    print WF "\t" . ${$hash{$key}}{$normal};
  }
  foreach my $tumor(@tumorSamples)
  {
    print WF "\t" . ${$hash{$key}}{$tumor};
  }
  print WF "\n";
}
close(WF);
##
####################################################################################################################################################################
exp <- read.table('mRNAmatrix.symbol.txt',header = T,sep = '\t',check.names = F)
colnames(exp)[1] <- 'id'
#
tcgaReplicateFilter = function(tsb, analyte_target=c("DNA","RNA"), decreasing=TRUE, analyte_position=20, plate=c(22,25), portion=c(18,19), filter_FFPE=FALSE, full_barcode=FALSE){
  # basically, user provide tsb and analyte_target is fine. If you
  # want to filter FFPE samples, please set filter_FFPE and full_barcode
  # all to TRUE, and tsb must have nchar of 28
  
  analyte_target = match.arg(analyte_target)
  # Strings in R are largely lexicographic
  # see ??base::Comparison
  
  # filter FFPE samples
  # provide by <http://gdac.broadinstitute.org/runs/sampleReports/latest/FPPP_FFPE_Cases.html> 
  if(full_barcode & filter_FFPE){
    ffpe = c("TCGA-44-2656-01B-06D-A271-08", "TCGA-44-2656-01B-06D-A273-01", 
             "TCGA-44-2656-01B-06D-A276-05", "TCGA-44-2656-01B-06D-A27C-26", 
             "TCGA-44-2656-01B-06R-A277-07", "TCGA-44-2662-01B-02D-A271-08", 
             "TCGA-44-2662-01B-02D-A273-01", "TCGA-44-2662-01B-02R-A277-07", 
             "TCGA-44-2665-01B-06D-A271-08", "TCGA-44-2665-01B-06D-A273-01", 
             "TCGA-44-2665-01B-06D-A276-05", "TCGA-44-2665-01B-06R-A277-07", 
             "TCGA-44-2666-01B-02D-A271-08", "TCGA-44-2666-01B-02D-A273-01", 
             "TCGA-44-2666-01B-02D-A276-05", "TCGA-44-2666-01B-02D-A27C-26", 
             "TCGA-44-2666-01B-02R-A277-07", "TCGA-44-2668-01B-02D-A271-08", 
             "TCGA-44-2668-01B-02D-A273-01", "TCGA-44-2668-01B-02D-A276-05", 
             "TCGA-44-2668-01B-02D-A27C-26", "TCGA-44-2668-01B-02R-A277-07", 
             "TCGA-44-3917-01B-02D-A271-08", "TCGA-44-3917-01B-02D-A273-01", 
             "TCGA-44-3917-01B-02D-A276-05", "TCGA-44-3917-01B-02D-A27C-26", 
             "TCGA-44-3917-01B-02R-A277-07", "TCGA-44-3918-01B-02D-A271-08", 
             "TCGA-44-3918-01B-02D-A273-01", "TCGA-44-3918-01B-02D-A276-05", 
             "TCGA-44-3918-01B-02D-A27C-26", "TCGA-44-3918-01B-02R-A277-07", 
             "TCGA-44-4112-01B-06D-A271-08", "TCGA-44-4112-01B-06D-A273-01", 
             "TCGA-44-4112-01B-06D-A276-05", "TCGA-44-4112-01B-06D-A27C-26", 
             "TCGA-44-4112-01B-06R-A277-07", "TCGA-44-5645-01B-04D-A271-08", 
             "TCGA-44-5645-01B-04D-A273-01", "TCGA-44-5645-01B-04D-A276-05", 
             "TCGA-44-5645-01B-04D-A27C-26", "TCGA-44-5645-01B-04R-A277-07", 
             "TCGA-44-6146-01B-04D-A271-08", "TCGA-44-6146-01B-04D-A273-01", 
             "TCGA-44-6146-01B-04D-A276-05", "TCGA-44-6146-01B-04D-A27C-26", 
             "TCGA-44-6146-01B-04R-A277-07", "TCGA-44-6146-01B-04R-A27D-13", 
             "TCGA-44-6147-01B-06D-A271-08", "TCGA-44-6147-01B-06D-A273-01", 
             "TCGA-44-6147-01B-06D-A276-05", "TCGA-44-6147-01B-06D-A27C-26", 
             "TCGA-44-6147-01B-06R-A277-07", "TCGA-44-6147-01B-06R-A27D-13", 
             "TCGA-44-6775-01C-02D-A271-08", "TCGA-44-6775-01C-02D-A273-01", 
             "TCGA-44-6775-01C-02D-A276-05", "TCGA-44-6775-01C-02D-A27C-26", 
             "TCGA-44-6775-01C-02R-A277-07", "TCGA-44-6775-01C-02R-A27D-13", 
             "TCGA-A6-2674-01B-04D-A270-10", "TCGA-A6-2674-01B-04R-A277-07", 
             "TCGA-A6-2677-01B-02D-A270-10", "TCGA-A6-2677-01B-02D-A274-01", 
             "TCGA-A6-2677-01B-02D-A27A-05", "TCGA-A6-2677-01B-02D-A27E-26", 
             "TCGA-A6-2677-01B-02R-A277-07", "TCGA-A6-2684-01C-08D-A270-10", 
             "TCGA-A6-2684-01C-08D-A274-01", "TCGA-A6-2684-01C-08D-A27A-05", 
             "TCGA-A6-2684-01C-08D-A27E-26", "TCGA-A6-2684-01C-08R-A277-07", 
             "TCGA-A6-3809-01B-04D-A270-10", "TCGA-A6-3809-01B-04D-A274-01", 
             "TCGA-A6-3809-01B-04D-A27A-05", "TCGA-A6-3809-01B-04D-A27E-26", 
             "TCGA-A6-3809-01B-04R-A277-07", "TCGA-A6-3810-01B-04D-A270-10", 
             "TCGA-A6-3810-01B-04D-A274-01", "TCGA-A6-3810-01B-04D-A27A-05", 
             "TCGA-A6-3810-01B-04D-A27E-26", "TCGA-A6-3810-01B-04R-A277-07", 
             "TCGA-A6-5656-01B-02D-A270-10", "TCGA-A6-5656-01B-02D-A274-01", 
             "TCGA-A6-5656-01B-02D-A27A-05", "TCGA-A6-5656-01B-02D-A27E-26", 
             "TCGA-A6-5656-01B-02R-A277-07", "TCGA-A6-5656-01B-02R-A27D-13", 
             "TCGA-A6-5659-01B-04D-A270-10", "TCGA-A6-5659-01B-04D-A274-01", 
             "TCGA-A6-5659-01B-04D-A27A-05", "TCGA-A6-5659-01B-04D-A27E-26", 
             "TCGA-A6-5659-01B-04R-A277-07", "TCGA-A6-6650-01B-02D-A270-10", 
             "TCGA-A6-6650-01B-02D-A274-01", "TCGA-A6-6650-01B-02D-A27A-05", 
             "TCGA-A6-6650-01B-02D-A27E-26", "TCGA-A6-6650-01B-02R-A277-07", 
             "TCGA-A6-6650-01B-02R-A27D-13", "TCGA-A6-6780-01B-04D-A270-10", 
             "TCGA-A6-6780-01B-04D-A274-01", "TCGA-A6-6780-01B-04D-A27A-05", 
             "TCGA-A6-6780-01B-04D-A27E-26", "TCGA-A6-6780-01B-04R-A277-07", 
             "TCGA-A6-6780-01B-04R-A27D-13", "TCGA-A6-6781-01B-06D-A270-10", 
             "TCGA-A6-6781-01B-06D-A274-01", "TCGA-A6-6781-01B-06D-A27A-05", 
             "TCGA-A6-6781-01B-06R-A277-07", "TCGA-A6-6781-01B-06R-A27D-13", 
             "TCGA-A7-A0DB-01C-02D-A272-09", "TCGA-A7-A0DB-01C-02R-A277-07", 
             "TCGA-A7-A0DB-01C-02R-A27D-13", "TCGA-A7-A13D-01B-04D-A272-09", 
             "TCGA-A7-A13D-01B-04R-A277-07", "TCGA-A7-A13D-01B-04R-A27D-13", 
             "TCGA-A7-A13E-01B-06D-A272-09", "TCGA-A7-A13E-01B-06R-A277-07", 
             "TCGA-A7-A13E-01B-06R-A27D-13", "TCGA-A7-A26E-01B-06D-A272-09", 
             "TCGA-A7-A26E-01B-06D-A275-01", "TCGA-A7-A26E-01B-06D-A27B-05", 
             "TCGA-A7-A26E-01B-06R-A277-07", "TCGA-A7-A26E-01B-06R-A27D-13", 
             "TCGA-A7-A26J-01B-02D-A272-09", "TCGA-A7-A26J-01B-02D-A275-01", 
             "TCGA-A7-A26J-01B-02D-A27B-05", "TCGA-A7-A26J-01B-02D-A27F-26", 
             "TCGA-A7-A26J-01B-02R-A277-07", "TCGA-A7-A26J-01B-02R-A27D-13", 
             "TCGA-B2-3923-01B-10D-A270-10", "TCGA-B2-3923-01B-10R-A277-07", 
             "TCGA-B2-3923-01B-10R-A27D-13", "TCGA-B2-3924-01B-03D-A270-10", 
             "TCGA-B2-3924-01B-03D-A274-01", "TCGA-B2-3924-01B-03D-A27A-05", 
             "TCGA-B2-3924-01B-03D-A27E-26", "TCGA-B2-3924-01B-03R-A277-07", 
             "TCGA-B2-3924-01B-03R-A27D-13", "TCGA-B2-5633-01B-04D-A270-10", 
             "TCGA-B2-5633-01B-04D-A274-01", "TCGA-B2-5633-01B-04D-A27A-05", 
             "TCGA-B2-5633-01B-04D-A27E-26", "TCGA-B2-5633-01B-04R-A277-07", 
             "TCGA-B2-5633-01B-04R-A27D-13", "TCGA-B2-5635-01B-04D-A270-10", 
             "TCGA-B2-5635-01B-04D-A274-01", "TCGA-B2-5635-01B-04D-A27A-05", 
             "TCGA-B2-5635-01B-04D-A27E-26", "TCGA-B2-5635-01B-04R-A277-07", 
             "TCGA-B2-5635-01B-04R-A27D-13", "TCGA-BK-A0CA-01B-02D-A272-09", 
             "TCGA-BK-A0CA-01B-02D-A275-01", "TCGA-BK-A0CA-01B-02D-A27B-05", 
             "TCGA-BK-A0CA-01B-02D-A27F-26", "TCGA-BK-A0CA-01B-02R-A277-07", 
             "TCGA-BK-A0CA-01B-02R-A27D-13", "TCGA-BK-A0CC-01B-04D-A272-09", 
             "TCGA-BK-A0CC-01B-04D-A275-01", "TCGA-BK-A0CC-01B-04D-A27B-05", 
             "TCGA-BK-A0CC-01B-04R-A277-07", "TCGA-BK-A0CC-01B-04R-A27D-13", 
             "TCGA-BK-A139-01C-08D-A272-09", "TCGA-BK-A139-01C-08D-A275-01", 
             "TCGA-BK-A139-01C-08D-A27B-05", "TCGA-BK-A139-01C-08D-A27F-26", 
             "TCGA-BK-A139-01C-08R-A277-07", "TCGA-BK-A139-01C-08R-A27D-13", 
             "TCGA-BK-A26L-01C-04D-A272-09", "TCGA-BK-A26L-01C-04D-A275-01", 
             "TCGA-BK-A26L-01C-04D-A27B-05", "TCGA-BK-A26L-01C-04D-A27F-26", 
             "TCGA-BK-A26L-01C-04R-A277-07", "TCGA-BK-A26L-01C-04R-A27D-13", 
             "TCGA-BL-A0C8-01B-04D-A271-08", "TCGA-BL-A0C8-01B-04D-A273-01", 
             "TCGA-BL-A0C8-01B-04D-A276-05", "TCGA-BL-A0C8-01B-04D-A27C-26", 
             "TCGA-BL-A0C8-01B-04R-A277-07", "TCGA-BL-A0C8-01B-04R-A27D-13", 
             "TCGA-BL-A13I-01B-04D-A271-08", "TCGA-BL-A13I-01B-04D-A276-05", 
             "TCGA-BL-A13I-01B-04R-A277-07", "TCGA-BL-A13I-01B-04R-A27D-13", 
             "TCGA-BL-A13J-01B-04D-A271-08", "TCGA-BL-A13J-01B-04D-A273-01", 
             "TCGA-BL-A13J-01B-04D-A276-05", "TCGA-BL-A13J-01B-04D-A27C-26", 
             "TCGA-BL-A13J-01B-04R-A277-07", "TCGA-BL-A13J-01B-04R-A27D-13")
    
    tsb = setdiff(tsb, tsb[which(tsb %in% ffpe)])
  }
  
  # find repeated samples
  sampleID = substr(tsb, start = 1, stop = 15)
  dp_samples = unique(sampleID[duplicated(sampleID)])
  
  if(length(dp_samples)==0){
    message("ooo Not find any duplicated barcodes, return original input..")
    tsb
  }else{
    uniq_tsb = tsb[! sampleID %in% dp_samples]
    dp_tsb = setdiff(tsb, uniq_tsb)
    
    add_tsb = c()
    
    # analyte = substr(dp_tsb, start = analyte_position, stop = analyte_position)
    # if analyte_target = "DNA"
    # analyte:  D > G,W,X
    if(analyte_target == "DNA"){
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "D") & !(all(analytes == "D"))){
          aliquot = mulaliquots[which(analytes == "D")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }else{
      # if analyte_target = "RNA"
      # analyte: H > R > T 
      for(x in dp_samples){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        analytes = substr(mulaliquots, 
                          start = analyte_position,
                          stop = analyte_position)
        if(any(analytes == "H") & !(all(analytes == "H"))){
          aliquot = mulaliquots[which(analytes == "H")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
          
        }else if(any(analytes == "R") & !(all(analytes == "R"))){
          aliquot = mulaliquots[which(analytes == "R")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }else if(any(analytes == "T") & !(all(analytes == "T"))){
          aliquot = mulaliquots[which(analytes == "T")]
          
          if (length(aliquot) != 1) {
            # Still have repeats
            # Remove the samples and add repeated id back to list
            dp_tsb = c(setdiff(dp_tsb, mulaliquots), aliquot)
          } else {
            # If have no repeats
            # Just remove samples from list and
            # add unique id to result list
            add_tsb = c(add_tsb, aliquot)
            dp_tsb = setdiff(dp_tsb, mulaliquots) 
          }
        }
        
      }
    }
    
    
    if(length(dp_tsb) == 0){
      message("ooo Filter barcodes successfully!")
      c(uniq_tsb, add_tsb)
    }else{
      # filter according to portion number
      sampleID_res = substr(dp_tsb, start=1, stop=15)
      dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
      
      for(x in dp_samples_res){
        mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
        portion_codes = substr(mulaliquots,
                               start = portion[1],
                               stop = portion[2])
        portion_keep = sort(portion_codes, decreasing = decreasing)[1]
        if(!all(portion_codes == portion_keep)){
          if(length(which(portion_codes == portion_keep)) == 1){
            add_tsb = c(add_tsb, mulaliquots[which(portion_codes == portion_keep)])
            dp_tsb = setdiff(dp_tsb, mulaliquots)
          }else{
            dp_tsb = setdiff(dp_tsb, mulaliquots[which(portion_codes != portion_keep)])
          }
          
        }
      }
      
      if(length(dp_tsb)==0){
        message("ooo Filter barcodes successfully!")
        c(uniq_tsb, add_tsb)
      }else{
        # filter according to plate number
        sampleID_res = substr(dp_tsb, start=1, stop=15)
        dp_samples_res = unique(sampleID_res[duplicated(sampleID_res)])
        for(x in dp_samples_res){
          mulaliquots = dp_tsb[substr(dp_tsb,1,15) == x]
          plate_codes = substr(mulaliquots,
                               start = plate[1],
                               stop = plate[2])
          plate_keep = sort(plate_codes, decreasing = decreasing)[1]
          add_tsb = c(add_tsb, mulaliquots[which(plate_codes == plate_keep)])
          dp_tsb = setdiff(dp_tsb, mulaliquots)
        }
        
        if(length(dp_tsb)==0){
          message("ooo Filter barcodes successfully!")
          c(uniq_tsb, add_tsb)
        }else{
          message("ooo barcodes ", dp_tsb, " failed in filter process, other barcodes will be returned.")
          c(uniq_tsb, add_tsb)
        }
      }
    }
  }
}
dup.filt <- tcgaReplicateFilter(colnames(exp),'RNA')
exp.filt <- exp[,dup.filt]
#
exp.filt1 <- aggregate(.~id,exp.filt,mean)
rownames(exp.filt1) <- exp.filt1$id
#exp.filt1 <- exp.filt1[,-1]
#
write.table(exp.filt1,'TCGA.exp.filt1.txt',sep = '\t',quote = F,row.names = F)
#
#exp.filt$ensemble <- unlist(str_split(exp.filt$id,"[.]",simplify=T))[,1]
#######################################################################################################################################################
#CIBERSORT
#
## ssGSEA
library(tidyr)
library(dplyr)
library(stringr)
#save results
setwd("/data/project-SunJL/yqxx0148/1-cibersort/mmeigui/normal/")
#CIBERSORT####
source("CIBERSORT.R")

result <- NULL
seedresult <- 0
res <- data.frame()
for (i in 1:5){
  set.seed(i)    ## seed=3,nsample=537
  results=CIBERSORT("LM22.txt", "TCGA.exp.filt1.txt", perm=100, QN=TRUE)
  filt <- read.table('CIBERSORT-Results.txt',header=T,sep="\t",check.names=F,row.names=1) 
  filt <- filt[filt$Pvalue <0.05,]
  res <- rbind(res,data.frame(seed=i,nsample = nrow(filt)))
  if (seedresult < nrow(filt)){
    seedresult <- nrow(filt)
    result <- read.table('CIBERSORT-Results.txt',header=T,sep="\t",check.names=F,row.names=1)
  }
}

write.table(result,'result.txt',sep = '\t',quote = F,row.names = F)

#

#
#barplot
input="CIBERSORT-Results-p0.05.txt" ###p<0.05
data <- read.table(input,header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf('infi_heatmap.pdf',height=12,width=18)
#png('infi_heatmap.png',height=800,width=1200)
par(las=1,mar=c(8,8,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F)
par(srt=90,xpd=T);text(a1,-0.04,colnames(data),adj=1,cex=1);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),
       col=col,pch=15,bty="n",cex=1.2)
dev.off()

#
input="CIBERSORT-Results.txt"
data <- read.table(input,header=T,sep="\t",check.names=F,row.names=1)
data=t(data)
col=rainbow(nrow(data),s=0.7,v=0.7)
pdf('CIBESORT.bar.pdf',height=10,width=18)
#png('infi_heatmap.png',height=800,width=1200)
par(las=1,mar=c(8,8,4,15))
a1 = barplot(data,col=col,yaxt="n",ylab="Relative Percent",xaxt="n")
a2=axis(2,tick=F,labels=F)
axis(2,a2,paste0(a2*100,"%"))
axis(1,a1,labels=F,tick = F)
par(xpd=T);text(115,-0.04,'samples of GSE20680 (n = 139)',adj=1,cex=1.7);par(srt=0)
ytick2 = cumsum(data[,ncol(data)])
ytick1 = c(0,ytick2[-length(ytick2)])
#text(par('usr')[2],(ytick1+ytick2)/2,rownames(data),cex=0.6,adj=0)
legend(par('usr')[2]*0.98,par('usr')[4],legend=rownames(data),
       col=col,pch=15,bty="n",cex=1.2)
dev.off()
########################################################################################################################################################
# 
#                            type     number                       （后面加的）label
# 1             T cells CD4 naive  0.0019100              T cells CD4 naive0.00191%
# 2           T cells gamma delta  0.0738263          T cells gamma delta0.0738263%
# 3                B cells memory  0.2145824               B cells memory0.2145824%
# 4    T cells regulatory (Tregs)  0.2637631   T cells regulatory (Tregs)0.2637631%
# 5  T cells CD4 memory activated  0.4219686 T cells CD4 memory activated0.4219686%
# 6     T cells follicular helper  0.4459848    T cells follicular helper0.4459848%
# 7          Mast cells activated  0.5056034         Mast cells activated0.5056034%
# 8                  Plasma cells  1.1309978                 Plasma cells1.1309978%
# 9                 B cells naive  1.3674974                B cells naive1.3674974%
# 10           NK cells activated  1.7255253           NK cells activated1.7255253%
# 11                  Eosinophils  1.7663008                  Eosinophils1.7663008%
# 12      Dendritic cells resting  2.2153714      Dendritic cells resting2.2153714%
# 13                  Neutrophils  2.6894470                   Neutrophils2.689447%
# 14    Dendritic cells activated  3.1519054    Dendritic cells activated3.1519054%
# 15               Macrophages M1  3.3892539               Macrophages M13.3892539%
# 16             NK cells resting  3.6298910              NK cells resting3.629891%
# 17                    Monocytes  5.5927492                    Monocytes5.5927492%
# 18                  T cells CD8  6.9256826                  T cells CD86.9256826%
# 19           Mast cells resting  9.1783956           Mast cells resting9.1783956%
# 20               Macrophages M0 17.6680276              Macrophages M017.6680276%
# 21   T cells CD4 memory resting 18.2295994  T cells CD4 memory resting18.2295994%
# 22               Macrophages M2 19.4117142              Macrophages M219.4117142%

# CIBERSORT-Results.txt中p>0.05
setwd("/data/project-SunJL/yqxx0148/1-cibersort/mmeigui/normal/rerun/")
library(ggplot2)
library(tidyverse)
input="normaltest.txt"        #
data <- read.table(input,sep="\t")
data <- t(data)

write.table(data,file = "data.txt",quote = F,sep = "\t")  #

data <- read.table("data-normal.txt",sep = "\t",header = T)
data$number <- round(data$number,3)

#lable
data <- data %>%
  mutate(
    label = case_when(
      rownames(data) <= 8 ~ paste0(type, data$number,"%"),
      rownames(data) <= 16 ~ paste0(type, data$number, "%\n"),
      T ~ paste0(type, "\n", number, "%")
    )
  )

data <- read.table("data.txt",sep = "\t",header = T)
data$number <- round(data$number,3)
write.table(data,file = "data-normal.txt",sep = "\t")
#
data$type<- factor(data$type,
                   levels = data$type)
##
p1 <- ggplot(data = data, aes(x = type, y = number, label = label)) +
  geom_col(aes(fill = type), width = 1, size = 0) +
  geom_col( aes(y = 0.3), fill = "white",width = 1,alpha = 0.2, size = 0 ) +
  geom_col( aes(y = 0), fill = "white", width = 1, alpha = 0.2, size = 0)
p1

##
p2 <-
  p1 +
  coord_polar() +
  theme_void() +
  scale_y_continuous(limits = c(-3, 20))   #
p2

##lable
p3 <-p2 +
  geom_text(data = . %>% filter(rownames(data) <= 22),nudge_y = 0.03,
            angle = 1,fontface = "bold",size = 1.5) +
  geom_text(data = . %>% filter( rownames(data)> 22),
            nudge_y = -0.05,
            color = "white",
            angle =1,
            fontface = "bold",
            size = 1.5)
p3

##
#####################################################################################################################################################
library(vioplot)                                                    #
rt = cibesort.final
# rt = cibesort
rt <- rt[,1:22] 
rt <- rt[,-5]
# loc <- which(apply(rt,2,function(x) sum(x==0) == length(x)))
# if(length(loc) > 0) {rt <- rt[,-loc]}
# a <- rt[order(as.numeric(substr(rownames(rt),9,10))),]
normal=58               # C1
case= 479          #C2 
# normal=58               # C1 
# case= 477

pdf("vioplot.v2.pdf",width = 14, height = 10)              #
#par(las=1,mar=c(13,6,3,3))
#tiff(file="vioplot.tiff",width = 15, height = 18,units = "in",res = 150)
x=c(1:ncol(rt))
y=c(1:ncol(rt))
plot(x,y,
     xlim=c(0,42),ylim=c(min(rt),max(rt)+0.04),
     main="",xlab="", ylab="Fraction",
     pch=21,
     col="white",
     xaxt="n")
TNpvalue <- list()
for(i in 1:ncol(rt)){
  normalData=rt[1:normal,i]
  tumorData=rt[(normal+1):(normal+case),i]
  vioplot(normalData,at=2*(i-1),lty=1,add = T,col = 'blue')
  vioplot(tumorData,at=2*(i-1)+1,lty=1,add = T,col = 'red')
  wilcoxTest=wilcox.test(normalData,tumorData)
  p=round(wilcoxTest$p.value,2)
  mx=max(c(normalData,tumorData),na.rm = T)
  lines(c(x=2*(i-1)+0.2,x=2*(i-1)+0.8),c(mx,mx))
  text(x=2*(i-1)+0.5, y=mx+0.04,labels=ifelse(p<0.05, paste0("p<0.05"), paste0("p=",p)), cex = 1)
  
  #text(seq(1,64,3),-0.03,labels=colnames(a),srt = 90)
  text(seq(1,43,2)[i],-0.03,xpd = NA,labels=colnames(rt)[i], font=1, cex = 1,srt = 75,pos=2)
  #print(ifelse(p<0.001, paste0("p<0.001"), paste0("p=",p)))
  TNpvalue[[i]] <- c(colnames(rt)[i],mean(tumorData),mean(normalData),p)
}
#axis(1,at=seq(1,63,3),labels=colnames(a),las = 2)
dev.off()

TNpvalue <- do.call(rbind,TNpvalue)
write.csv(TNpvalue,file="TNpvalue.csv")
###################################################################################################################################################
## wgcna
## TCGA.exp.filt1.txt
##
setwd("/data/project-SunJL/yqS095/")
library(tidyr)
library(dplyr)
library(stringr)
#
exp.final <- read.table('TCGA.exp.filt1.txt',sep = '\t',check.names = F,header = T,row.names=1)
#
sample <- read.table('sample.txt',sep='\t',check.names=F,header=T)
group_list <- factor(sample$group,levels = c("normal","tumor"))
exprSet <- exp.final
#log2
#
exprSet <- apply(exp.final,2,function(x){log2(x+1)})
#
pvalue <-0.05
logFoldChange <- 1
dat <- exprSet
design = model.matrix(~0+group_list, data=group_list)
colnames(design) <- levels(group_list)
rownames(design) <- colnames(dat)
design
#
contrast.matrix <- makeContrasts(tumor-normal, levels = design)#Case比Control
fit <- lmFit(dat, design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit2 <- eBayes(fit2)
options(digits = 4)#
allDiff=topTable(fit2,coef=1,number=Inf)
allDiff <- na.omit(allDiff)
write.table(cbind(Symbol=rownames(allDiff),allDiff),file="limmaOut.txt",sep="\t",row.names = F,quote = F)

##
diffSig = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange | allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffSig),diffSig), file="diffSig.txt",sep="\t",row.names = F,quote = F)
diffUp = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC>logFoldChange)),]
write.table(cbind(Symbol=rownames(diffUp),diffUp), file="up.txt",sep="\t",row.names = F,quote = F)
diffDown = allDiff[(allDiff$P.Value < pvalue & (allDiff$logFC<(-logFoldChange))),]
write.table(cbind(Symbol=rownames(diffDown),diffDown), file="down.txt",sep="\t",row.names = F,quote = F)
##
library(ggplot2)
pvalue <-0.05
logFC <- 1
allDiff <- read.table('limmaOut.txt',header = T,check.names = F,row.names = 1,sep = '\t')
allDiff$Significant <- ifelse(allDiff$P.Value<pvalue & abs(allDiff$logFC)>= logFC,
                              ifelse(allDiff$logFC> logFC,'up','down'),'no')
mycol <- c("#3CB371","#3D3D3D","#FF4500")
pdf(file="Fig1.MDD.DEGs.pdf",width=6,height=5.5)
p <- ggplot(allDiff, aes(logFC, -log10(P.Value), colour= Significant))+
  geom_point(size=1.2,alpha=0.4)+theme_bw()+
  scale_color_manual(values = mycol,name='Significant')+
  labs(title="DEGs",x="Log2FC",y="-log10 (P.value)")+
  geom_hline(yintercept = -log10(pvalue),linetype=3,lwd = 1)+
  geom_vline(xintercept = c(-logFC, logFC), linetype=3,lwd = 1)+
  theme(axis.text.x = element_text(size = 15),axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),text = element_text(size = 15),
        plot.title = element_text(hjust = 0.5))
print(p)
dev.off()
##
library(pheatmap)
diff <- rownames(diffSig)
diffexp <- dat[diff,]
annotation_col <- data.frame(Type = factor(group_list,levels = c("normal","tumor")))
rownames(annotation_col) <- colnames(diffexp)
color.key <- c("#3300CC","#3399FF","white","#FF3333","#CC0000")
pdf(file = 'Fig2.DEGs.MDD.heatmap.pdf',width = 12,height = 10,onefile=FALSE)  
pheatmap(diffexp,cellwidth = 0.5,cellheight = 0.2,fontsize_row=7,
         method="spearman", 
         scale='row',
         cluster_rows=TRUE,
         cluster_cols=FALSE,
         color = colorRampPalette(color.key)(50),
         show_colnames=FALSE,show_rownames =FALSE,
         annotation_col = annotation_col,
         #treeheight_row = "0",treeheight_col = "0",
         border_color = "NA")
dev.off()
##
diffsig1 <- diffSig
diffsig1 <- data.frame(names = row.names(diffsig1),diffsig1) 
diffsig2 <- diffsig1[1]       
rownames(diffsig2) <- NULL    
names(diffsig2) <- c("V1")    
write.table(diffSig2,file = "diffSigname.txt",quote=FALSE,sep="\t")
#
exp.diffsigname <- exp.final1[diffsigname$V1,]
write.table(exp.diffsigname,file = "exp.diffsignames.txt",sep = "\t",quote = F)
######################################################################################################################################################
setwd("/data/project-SunJL/yqxx0148/0-rerun/3wgcna/")
cibesort.final <- read.table('CIBERSORT-Results.txt',header=T,sep="\t",check.names=F,row.names=1) 
#cibesort <- read.table('CIBERSORT-Results.txt',sep = '\t',
#                       check.names = F,header = T,row.names = 1)
#cibesort.final <- rbind(cibesort[substr(rownames(cibesort),14,16)=='11A',],
#                        cibesort[substr(rownames(cibesort),14,16)=='01A',])
#cibesort.final 
exp <- read.table('TCGA.exp.filt1.txt',sep = '\t',check.names = F,header = T)
#tmp <- cbind(exp$id,exp[,substr(colnames(exp),14,15)=='11'],
#             exp[,substr(colnames(exp),14,16)=='01A'])
#write.table(tmp,'TCGA.exp.filt2.txt',sep = '\t',quote = F,row.names = F)
exp.filt <- exp[,c('id',rownames(cibesort.final))]
#exp.filt
rownames(exp.filt) <- exp.filt$id

#tumor <- rownames(cibesort[substr(rownames(cibesort),14,16)=='01A',])
tumor <- rownames(cibesort.final[59:537,])
#cibesort.final

#WGCNA
dataExpr <-data.frame(t(exp.filt[,tumor]),check.names = F)
dim(dataExpr)
library(WGCNA)
#
gsg = goodSamplesGenes(dataExpr, verbose = 3);
if (!gsg$allOK){
  # Optionally, print the gene and sample names that were removed:
  if (sum(!gsg$goodGenes)>0) 
    printFlush(paste("Removing genes:", 
                     paste(names(dataExpr)[!gsg$goodGenes], collapse = ",")));
  if (sum(!gsg$goodSamples)>0) 
    printFlush(paste("Removing samples:", 
                     paste(rownames(dataExpr)[!gsg$goodSamples], collapse = ",")));
  # Remove the offending genes and samples from the data:
  dataExpr = dataExpr[gsg$goodSamples, gsg$goodGenes]
}
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
dim(dataExpr)
save(dataExpr,file = 'dataExpr.Rdata')
sampleTree = hclust(dist(dataExpr), method = "complete")
pdf(file = "Fig1.Sample_cluster.pdf", width = 7, height = 5)
par(mar = c(0,4,2,0),xpd=F)
plot(sampleTree, main = "Sample clustering to detect outliers", cex=0.1,sub="", xlab="",cex.main=1.5)
#abline(h = 91, col = "red")
dev.off()
#
dataTraits <- cibesort.final[tumor,-c(4,5,12,20,23,24,25)]
head(dataTraits)
sampleTree2 = hclust(dist(dataExpr), method = "complete")
traitColors = numbers2colors(dataTraits, signed = FALSE)
pdf("Fig2.cluster.dendrogram.pdf",width=10,height=8)
plotDendroAndColors(sampleTree2, traitColors,
                    groupLabels = colnames(dataTraits),cex.dendroLabels=0.1,
                    main = "Sample dendrogram and trait heatmap")
dev.off()

###power
enableWGCNAThreads()   
powers =seq(from = 1, to=20, by=1)  #1:20
sft = pickSoftThreshold(dataExpr, powerVector = powers, verbose = 5)
save(sft,file='sft.Rdata')
pdf('Fig3.Soft.Threshold.pdf',width = 10,height = 6)
par(mfrow = c(1,2))
cex1 = 0.9
###
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=cex1,col="red") 
###
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")
dev.off()

###
softPower =sft$powerEstimate #power值
#
adjacency = adjacency(dataExpr, power = softPower)
softPower
###TOM
TOM = TOMsimilarity(adjacency)
dissTOM = 1-TOM
save(dissTOM,file='dissTOM.Rdata')
#load('dissTOM.Rdata')
###
geneTree = hclust(as.dist(dissTOM), method = "average")
save(geneTree,file='geneTree.Rdata')

###
minModuleSize =  100  #
dynamicMods = cutreeDynamic(dendro = geneTree, distM = dissTOM,
                            deepSplit = 2, pamRespectsDendro = FALSE,
                            minClusterSize = minModuleSize);
save(dynamicMods,file='dynamicMods.Rdata')
# load('dynamicMods.Rdata')
table(dynamicMods)
dynamicColors = labels2colors(dynamicMods)
table(dynamicColors)
pdf(file="Fig4.DynamicTree.pdf",width=6.5,height=4.5)
plotDendroAndColors(geneTree, dynamicColors, "Dynamic Tree Cut",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

###
MEList = moduleEigengenes(dataExpr, colors = dynamicColors)
MEs = MEList$eigengenes
#
moduleColors = MEList$validColors
table(moduleColors)
# MEDiss = 1-cor(MEs);
# METree = hclust(as.dist(MEDiss), method = "average")
# pdf('Fig5.merge.cluster.pdf',width = 9,height = 5)
# plot(METree, main = "Clustering of module eigengenes",
#      xlab = "", sub = "")
# MEDissThres = 0.5 
# abline(h=MEDissThres, col = "red")
# dev.off()
# 
# ###
# merge = mergeCloseModules(dataExpr, dynamicColors, cutHeight = MEDissThres, verbose = 3)
# mergedColors = merge$colors
# mergedMEs = merge$newMEs
# pdf(file="Fig6.DynamicTree.pdf",width=6.5,height=5)
# plotDendroAndColors(geneTree, cbind(dynamicColors,mergedColors),c("Dynamic Tree Cut","Merged dynamic"),
#                     dendroLabels = FALSE, hang = 0.03,
#                     addGuide = TRUE, guideHang = 0.05,
#                     main = "Gene dendrogram and module colors")
# dev.off()
# moduleColors = mergedColors
# save(moduleColors,file = 'moduleColors.Rdata')
# table(moduleColors)
# colorOrder = c("grey", standardColors(30))
# moduleLabels = match(moduleColors, colorOrder)-1
# MEs = mergedMEs

###
nGenes = ncol(dataExpr)
nSamples = nrow(dataExpr)
moduleTraitCor = cor(MEs, dataTraits, use = "p")
moduleTraitPvalue = corPvalueStudent(moduleTraitCor, nSamples)
textMatrix = paste(signif(moduleTraitCor, 2), "\n(",
                   signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) = dim(moduleTraitCor)
write.table(moduleTraitPvalue,'moduleTraitPvalue.txt',sep = '\t',quote = F,row.names = T)
write.table(moduleTraitCor,'moduleTraitCor.txt',sep = '\t',quote = F,row.names = T)
pdf(file="Fig7.Module_trait.pdf",width=10,height=13)
par(mar = c(9, 10, 3, 3))


#
labeledHeatmap(Matrix = moduleTraitCor,
               xLabels = colnames(dataTraits),
               yLabels = names(MEs),
               ySymbols = names(MEs),
               colorLabels = FALSE,
               colors = blueWhiteRed(500),
               textMatrix = textMatrix,
               setStdMargins = FALSE,
               cex.text =0.7,
               main = paste("Module-trait relationships"))
dev.off()

###MM,GS
#modNames = substring(names(MEs), 3)
#modNames <- c('brown','greenyellow','yellow','cyan','grey60','pink')
modNames <- unique(moduleColors)
geneModuleMembership = as.data.frame(cor(dataExpr, MEs, use = "p"))
MMPvalue = as.data.frame(corPvalueStudent(as.matrix(geneModuleMembership), nSamples))
#names(geneModuleMembership) = paste("MM", modNames, sep="")
#names(MMPvalue) = paste("p.MM", modNames, sep="")
traitNames=colnames(dataTraits)
geneTraitSignificance = as.data.frame(cor(dataExpr, dataTraits, use = "p"))
GSPvalue = as.data.frame(corPvalueStudent(as.matrix(geneTraitSignificance), nSamples))
names(geneTraitSignificance) = paste("GS.", traitNames, sep="")
names(GSPvalue) = paste("p.GS.", traitNames, sep="")

#
pdf('Fig8.Gene.trait.blue.pdf',width = 15,height = 22)
par(mfrow=c(6,3),mar=c(6,5,5,3))
for (i in colnames(geneTraitSignificance)){
  moduleGenes = moduleColors=='blue'
  verboseScatterplot(abs(geneModuleMembership[moduleGenes, 'MEblue']),
                     abs(geneTraitSignificance[moduleGenes,i]),
                     xlab = ("Module Membership in blue module"),
                     ylab = paste("Gene significance"),
                     main = paste0('Trait ',substr(i,4,nchar(i)),"\n"),
                     cex.main = 1.5, cex.lab = 1.5, cex.axis = 1.5, col = 'blue')
  abline(v=0.6,h=0.2,col="red")
}
dev.off()

###
for (mod in 1:nrow(table(moduleColors)))
{  
  modules = names(table(moduleColors))[mod]
  probes = colnames(dataExpr)
  inModule = (moduleColors == modules)
  modGenes = probes[inModule]
  write.table(modGenes, file =paste0('WGCNA.',modules,".txt"),sep=",",row.names=F,col.names=F,quote=F)
}

###GS_MM
probes = colnames(dataExpr)
geneInfo0 = data.frame(probes= probes,
                       moduleColor = moduleColors)
for (Tra in 1:ncol(geneTraitSignificance))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneTraitSignificance[,Tra],
                         GSPvalue[, Tra])
  names(geneInfo0) = c(oldNames,names(geneTraitSignificance)[Tra],
                       names(GSPvalue)[Tra])
}

for (mod in 1:ncol(geneModuleMembership))
{
  oldNames = names(geneInfo0)
  geneInfo0 = data.frame(geneInfo0, geneModuleMembership[,mod],
                         MMPvalue[, mod])
  names(geneInfo0) = c(oldNames,names(geneModuleMembership)[mod],
                       names(MMPvalue)[mod])
}
geneOrder =order(geneInfo0$moduleColor)
geneInfo = geneInfo0[geneOrder, ]
write.table(geneInfo, file = "GS_MM.csv",sep=",",row.names=F)
######################################################################################################################################################
##############################################################################-----lasso-svm-rf-lr-xgb-gbdt-----###################################################################################
AUC 4 ，seed = 18
#
library(glmnet)
library(tidyverse)
source('msvmRFE.R')
library(e1071)
library(caret)
library(randomForest)
library(RColorBrewer)
library(VennDiagram)
library(sigFeature)
library(e1071)
library(dplyr)
library(pROC)
Sys.setenv(LANGUAGE = "en") #
options(stringsAsFactors = FALSE) #chr,factor
rtt <- read.table("dmrgmachine.txt",row.names = 1,header = T,check.names = F,stringsAsFactors = F) 
sam <- rownames(rtt) 
rtt$group <- factor(ifelse(rtt$group=='no',0,1))
mycol <- RColorBrewer::brewer.pal(5,'Set1')
createFolds <- function(strat_id, k, seed = 123456) {
  set.seed(seed)
  if(k > length(strat_id)) {
    k <- length(strat_id)
  } 
  perm <- sample(length(strat_id))
  strat_id <- factor(strat_id, levels=sample(unique(strat_id)))
  strat_order <- order(strat_id[perm])
  num_item_per_fold_ceil <- ceiling(length(strat_id) / k)
  fold_ids <- rep(seq_len(k), times= num_item_per_fold_ceil)
  fold_ids <- fold_ids[seq_along(strat_id)]
  folds <- split(perm[strat_order], fold_ids)
  names(folds) <- paste0("Fold", seq_len(k))    
  return(folds)
}
#
n.fold <- 5 
seed = 12345678 
fold <- createFolds(sam,n.fold,seed = seed)
# 
roc <- data.frame()
letters1 <- letters
i <- 1
while(i<100){
  letters[i] <- paste(sample(letters1,2,replace = FALSE),collapse = '')
  i <- i+1
}
letters <- unique(letters)
set.seed(1)
letters <- sample(letters,44,replace = FALSE)

#####
#GBDT——https://zhuanlan.zhihu.com/p/25805870
#LASSO
#ROC
roc <- data.frame()
for (i in 1:n.fold){
  test_pred_RF <- pred_svm <- NULL
  train_sam <- sam[-fold[[i]]] 
  test_sam <- setdiff(sam,train_sam) 
  # 
  train_dat <- rtt[train_sam,]
  # 
  test_dat <- rtt[test_sam,]
  set.seed(seed)             
  x <- as.matrix(train_dat[,setdiff(colnames(train_dat),"group")])
  y <- train_dat$group
  cvfit = cv.glmnet(x, y, 
                    nfold=10, #10-fold cross-validation
                    family = "binomial", group.measure = "class")
  myCoefs <- coef(cvfit, s="lambda.min")
  lasso_fea <- rownames(coef(cvfit, s = 'lambda.min'))[coef(cvfit, s = 'lambda.min')[,1]!= 0] 
  lasso_fea <- setdiff(lasso_fea,"(Intercept)")
  if (length(lasso_fea)<=1) {
    next
  }
  
  #lasso features RF
  cat("RF",i,"...\n")
  set.seed(seed)   
  RFdata <- train_dat[,c("group",lasso_fea)]
  colnames(RFdata) <- gsub('-','_', colnames(RFdata))
  RF_model <- randomForest(group ~ ., 
                           data = RFdata,
                           ntree = 1000, 
                           nPerm = 50, 
                           mtry = floor(sqrt(ncol(RFdata)-1)), 
                           proximity = T,
                           importance1 = T)
  
  # 
  pred_RF <-data.frame(prob = predict(RF_model, newdata = test_dat,group="response"),
                       group = test_dat$group,
                       stringsAsFactors = F)
  cat("\n")
  
  ## svm
  set.seed(seed)    ##seed
  svm_model = svm(group ~ ., 
                  data=train_dat,
                  method = "svmRadial",
                  trControl = train_dat,
                  tuneLength = 8,
                  metric = "ROC")
  pred_svm <- data.frame(prob = predict(svm_model, newdata = test_dat,type="response"),
                         group = test_dat$group,
                         stringsAsFactors = F)
  pdf("llrxgs.pdf")
  #ROC
  library(ggplot2)
  library(pROC)
  pred_RF_prob <- pred_RF$prob
  pred_RF_prob <- as.numeric(pred_RF_prob)
  pred_RF_prob <- pred_RF_prob-1
  RF.roc <- plot.roc(pred_RF$group,pred_RF_prob,ylim=c(0,1),xlim=c(1,0),
                     smooth=F, 
                     ci=TRUE, 
                     main="",
                     col=mycol[2],
                     lwd=2, 
                     legacy.axes=T,
                     print.auc = F,
                     add = T)
  pred_svm_prob <- pred_svm$prob
  pred_svm_prob <- as.numeric(pred_svm_prob)
  pred_svm_prob <- pred_svm_prob-1
  svm.ROC <- plot.roc(pred_svm$group, pred_svm_prob,ylim=c(0,1),xlim=c(1,0),
                      smooth=F, 
                      ci=TRUE, 
                      main="",
                      col=mycol[5],
                      lwd=2, 
                      legacy.axes=T,
                      print.auc = F,
                      add = T)
  roc <- rbind(roc,data.frame(nfold = i,seed=seed,roc1=LR.roc$auc,roc2=RF.roc$auc,
                              roc3=XGB.ROC$auc,roc4=gbm.ROC$auc,roc5=svm.ROC$auc))
  legend.label <- c("AUC",paste0("RF: ",round(RF.roc$auc,3)),paste0("svm: ",round(svm.ROC$auc,3)))
  legend("bottomright", 
         legend = legend.label,
         col = mycol[1:3],
         lwd = 2,
         bty="n")
  invisible(dev.off())
}
write.table(roc,file = "roc-rfsvm.txt", sep = "\t", quote = FALSE)

##Accuracy rate
# pred_LR <- data.frame(prob = predict(LR_model, newdata = test_dat,group="response"),
# group = test_dat$group,
# stringsAsFactors = F)
#> pred_LR
#                  prob group
# GSM1784990  178.57264     0
# GSM1785006 -197.33614     0
# GSM1785013  -84.54876     0
# GSM1785021 -152.73738     0
# GSM1785029  -26.38379     0
# GSM1785050 -172.94766     0
# GSM1785054 -144.45260     0
# GSM1784988  186.13555     1
# GSM1784995  -43.57132     1
# GSM1784997  -49.25569     1
# GSM1785023   83.09465     1
# GSM1785045  130.40823     1

#lasso
library(glmnet)
fit = glmnet(x, y, family = "binomial", alpha = 1, lambda = NULL)
pdf("A_lasso.pdf", width = 5, height = 8)
plot(fit, xvar = "dev", label = TRUE)
dev.off()
cvfit = cv.glmnet(x, y, 
                  nfold=10, 
                  family = "binomial", type.measure = "class")
plot(cvfit)
cvfit$lambda.min
myCoefs <- coef(cvfit, s="lambda.min");
lasso_fea <- myCoefs@Dimnames[[1]][which(myCoefs != 0 )]
(lasso_fea <- lasso_fea[-1])
write.csv(lasso_fea,"feature_lasso.csv")

###############################################################################################################################################################
redmrgs <- read.table("top10.txt")
redmrgs$group <- factor(ifelse(redmrgs$group=='no',0,1))
redmrgs1 <- redmrgs
hub <- rtt
hubexp <- redmrgs1
rocdata <- data.frame(Sample = rownames(hubexp),
                      exp = hubexp$OGDH,
                      Type = hubexp$group) 
#mycol <- brewer.pal(10,'Set3')
mycol <- brewer.pal(10,'Set2')
pdf('Fig19.markerROC.pdf',width = 7,height = 7)
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, 
              main="",
              #print.thres="best", 
              col=mycol[1],
              lwd=2, 
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(hubexp)[1],' AUC : ',format(as.numeric(x$auc),digits=3))
for (i in colnames(hubexp[,2:ncol(hubexp)])){
  j=j+1
  rocdata <- data.frame(Sample = rownames(hubexp),
                        exp = hubexp[,i],
                        Type = hubexp$group) 
  x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
                smooth=F, 
                main="",
                #print.thres="best", sensitivity+ specificity
                col=mycol[j],
                lwd=2, 
                legacy.axes=T,add=T)
  auc.test <- c(auc.test,
                paste0(i,' AUC : ',format(as.numeric(x$auc),digits=3)))
}
legend(0.6,0.2, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()
auc.test  
write.table(auc.test,file = "auctest.txt",quote = F,sep = "\t")

# 
hub <- rtt
hubexp <- redmrgs1
rocdata <- data.frame(Sample = rownames(hubexp),
                      exp = hubexp$IMPA2,             
                      Type = hubexp$group) 
#mycol <- brewer.pal(10,'Set3')
mycol <- brewer.pal(10,'Set2')
pdf('IMPA2-ROC.pdf',width = 7,height = 7)           
x <- plot.roc(rocdata[,3],rocdata[,2],ylim=c(0,1),xlim=c(1,0),
              smooth=F, 
              main="",
              #print.thres="best", ensitivity+ specificity
              col=mycol[1],
              lwd=2, 
              grid=c(0.1,0.1),
              grid.col=c("black", "black"),
              legacy.axes=T)
j=1
auc.test <- paste0(colnames(hubexp)[8],' AUC : ',format(as.numeric(x$auc),digits=3)) 
legend(0.6,0.2, auc.test,lwd=2,bty="n",col=mycol,cex = 1.3)
dev.off()
####################################################################################################################################################
setwd("/data/project-SunJL/yq159-6/geo-gse70493/test/1/10-roc-buchong/")
rm(list=ls())
library("pROC")
library(ggplot2)
library(ggthemes)
library(verification)

data<-read.table("10-ROC.txt",fileEncoding="UTF-8",header=T,sep="\t")               
data<-na.omit(data)
pdf(file="OGDH_ROC.pdf",width=12,height=9,,family="GB1",onefile=FALSE)
roc.list <-roc(group~OGDH,plot=TRUE,ci=T,print.auc=TRUE,levels=c("0","1"),data = data)
levels(data$group)<-c('0','1')
roc_area<-roc.area(as.numeric(as.vector(data$group)), roc.list$predictor)
g <- ggroc(roc.list,legacy.axes = TRUE,colour = "red",size = 2)
g<-g + theme_base()+ xlab("1-Specificity") + ylab("Sensitivity")+ggtitle("OGDH_ROC curve")+theme(plot.margin=unit(rep(3,4),'lines'),plot.title = element_text(size=22,face="bold",hjust = 0.5),axis.title=element_text(size=20, face="bold"),axis.text=element_text(size=18, face="bold")) +geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype="dashed",colour="blue")+geom_text(aes(x=0.7,y=0.4,label=paste("AUC:",round(roc.list$auc,3),sep="")),size =10)+geom_text(aes(x=0.7,y=0.35,label=paste("P<",round(roc_area$p.value,8),sep="")),size =10)
print(g)
dev.off()

##
rm(list=ls())
library("pROC")
library(ggplot2)
library(ggthemes)
library(verification)

data<-read.table("10-ROC.txt",fileEncoding="UTF-8",header=T,sep="\t")
data<-na.omit(data)
pdf(file="Logistic_Regression_model_ROC.pdf",width=12,height=9,,family="GB1",onefile=FALSE)
LR_model <-glm(as.factor(group)~OGDH+PLA2G4A+PLA2G15+PTGIS+PLA2G2C+UGT2B15+BAAT+IMPA2+FMO5+MTHFD1L,family=binomial(link = "logit"),data = data)
glm.probs <- predict(LR_model,data[,c("OGDH","PLA2G4A","PLA2G15","PTGIS","PLA2G2C","UGT2B15","BAAT","IMPA2","FMO5","MTHFD1L")],type="response")
data$model_score<-glm.probs
roc.list <-roc(group~model_score,plot=TRUE,ci=T,print.auc=TRUE,levels=c("0","1"),data = data)
levels(data$group)<-c('0','1')
roc_area<-roc.area(as.numeric(as.vector(data$group)), roc.list$predictor)
g <- ggroc(roc.list,legacy.axes = TRUE,colour = "green",size = 2)
g<-g + theme_base()+ xlab("1-Specificity") + ylab("Sensitivity")+ggtitle("ROC curve for Logistic Regression model")+theme(plot.margin=unit(rep(3,4),'lines'),plot.title = element_text(size=22,face="bold",hjust = 0.5),axis.title=element_text(size=20, face="bold"),axis.text=element_text(size=18, face="bold")) +geom_segment(aes(x = 0, xend = 1, y = 0, yend = 1), linetype="dashed",colour="blue")+geom_text(aes(x=0.7,y=0.4,label=paste("AUC:",round(roc.list$auc,3),sep="")),size =10)+geom_text(aes(x=0.7,y=0.35,label=paste("P<",round(roc_area$p.value,10),sep="")),size =10)
print(g)
dev.off()

#####################################################################################################################################################################
pcatest <- read.table("top10.txt")
ord <- colnames(pcatest)
ord <- ord[-1]

# expr_df <- read.table(file='top10pca.txt',stringsAsFactors = FALSE)
# colnames(expr_df) <- expr_df[1,]
# rownames(expr_df) <- expr_df[,1]
# expr_df <- expr_df[-1,]
# expr_df1 <- expr_df[,-1]
# 
# expr_df1[1:3,1:4]

meta_df <- read.table(file='top10group.txt',stringsAsFactors = FALSE)
colnames(meta_df) <- meta_df[1,]
meta_df <- meta_df[-1,]
#
head(meta_df, n=3)
#
pca.results <- prcomp(pcatest[,-1], center = TRUE, scale. = FALSE)

#
mycol <- c("#223D6C","#D20A13","#088247","#FFD121","#11AA4D","#58CDD9","#7A142C","#5D90BA","#431A3D","#91612D","#6E568C","#E0367A","#D8D155","#64495D","#7CC767")

#
#install.packages("devtools")
#library(devtools)
#devtools::install_github("GuangchuangYu/yyplot")
#devtools::install_github('fawda123/ggord')
library(ggplot2)
library(plyr)
library(ggord)
#library(yyplot)

source('./geom_ord_ellipse.R') 
#ggord PCA
ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
      ellipse = FALSE, 
      size = 2, 
      alpha=0.5, 
      #
      cols = mycol[1:length(unique(meta_df$group))],
      arrow = NULL,txt = NULL) + 
  theme(panel.grid =element_blank()) + 
  
  #yyplot
  geom_ord_ellipse(ellipse_pro = .95, 
                   size=1.5, 
                   lty=1 ) 
#pdf
ggsave("PCA_classic.pdf", width = 6, height = 6)

#
ggord(pca.results, grp_in = meta_df$group, repel=TRUE,
      alpha = 0.6,
      #alpha_el = 0.3,
      ellipse_pro = 0.95,
      size = 2,
      #cols = mycol[1:length(unique(meta_df$group))],
      arrow=0.2, 
      vec_ext = 5,
      veccol="brown",
      txt=3) + 
  theme(panel.grid =element_blank()) + 
  
  geom_ord_ellipse(ellipse_pro = .95, 
                   color='darkgrey', 
                   size=0.5, 
                   lty=2 ) + 
  geom_ord_ellipse(ellipse_pro = .98, 
                   #color='grey', 
                   size=0.5, lty=2 ) 
#pdf
ggsave("PCA_arrow.pdf", width = 6, height = 6)
##########################################################################################################
rm(list = ls())
save.fig=F
suppressMessages(library(igraph))
suppressMessages(library(Seurat))
suppressMessages(library(SingleR))
suppressMessages(library(Matrix))
suppressMessages(library(dplyr))
suppressMessages(library(sctransform))
suppressMessages(library(scater))
suppressMessages(library(stringr))
suppressMessages(library(progress))
suppressMessages(library(edgeR))
suppressMessages(library(GGally))
suppressMessages(library(clustree))
suppressMessages(library(tidytree))
suppressMessages(library(kableExtra))
suppressMessages(library(gdata))
suppressMessages(library(enrichplot))
suppressMessages(library(clusterProfiler))
suppressMessages(library(UpSetR))
suppressMessages(library(pheatmap))
# suppressMessages(library(YQSYtools))
suppressMessages(library(openxlsx))
suppressMessages(library(DropletUtils))
suppressMessages(library(scDblFinder))

dir.create('singlecell')
dir.create('singlecell/0.rawdata')

dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
dirlist <- dir(path = 'singlecell/0.rawdata',full.names = T)
path='singlecell'
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dirlist[which(samplename %in% clinical$X.Sample_geo_accession)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]


object.list <- list()
pb <- progress_bar$new(format = "  over [:bar] :percent time :elapsed", 
                       total = length(dirlist), clear = FALSE, width = 60)

for (i in 1:length(samplename)) {
  file = dirlist[i]
  # 
  if (!file.exists(paste0(path, "/1.QC/"))) {
    dir.create(paste0(path, "/1.QC/"))
  }
  if (!file.exists(paste0(path, "/1.QC/", samplename[i]))) {
    dir.create(paste0(path, "/1.QC/", samplename[i]))
  }
  message("")
  message("read file: ", samplename[i])
  # read data-read expression matrix-BD expression matrix csv file
  data <- t(read.csv(file,header = T,row.names = 1))
  
  object <- CreateSeuratObject(counts = data, project = samplename[i], min.cells = 3,         min.features = 200)
  
  # add the column of mtochondria percentage features
  # PercentageFeatureSet
  object[["percent.mt"]] <- PercentageFeatureSet(object, pattern = "^MT-")
  # QC-statistics of mt.percent per cell 
  MTloc <- grep('^MT-',rownames(object))
  print(summary(object[["percent.mt"]]))
  mt.percnet_summary <- summary(object[["percent.mt"]])
  capture.output(mt.percnet_summary, file = paste0(path, "/1.QC/", samplename[i], 
                                                   "/mt.percent_statitics.txt"))
  
  # add the column of HB percentage features 
  HB.genes_total <- c("HBA1", "HBA2", "HBB", "HBD", "HBG1", "HBQ1")
  HB_m <- match(HB.genes_total, rownames(object@assays$RNA))
  HB.genes <- rownames(object@assays$RNA)[HB_m]
  HB.genes <- HB.genes[!is.na(HB.genes)]
  object[["percent.HB"]] <- PercentageFeatureSet(object, features = HB.genes)
  # QC-statistics of mt.percent per cell
  print(summary(object[["percent.HB"]]))
  HB.percnet_summary <- summary(object[["percent.HB"]])
  capture.output(HB.percnet_summary, file = paste0(path, "/1.QC/", samplename[i], 
                                                   "/HB.percent_statitics.txt"))
  
  # QC-statistics of UMI counts per cell count
  summary(colSums(object))
  UMI_counts_summary <- summary(colSums(object))
  capture.output(UMI_counts_summary, file = paste0(path, "/1.QC/", samplename[i], 
                                                   "/UMI_counts_statitics.txt"))
  
  p <- hist(colSums(object), breaks = 100)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts.pdf"))
  # QC-histogram of UMI counts distribution
  hist(colSums(object), breaks = 100, main = "Distribution of UMI counts", xlab = "UMI counts", 
       ylab = "Cell number", col = "blue", xlim = c(0, max(p$breaks) * 1.1), ylim = c(0, 
                                                                                      max(p$counts) * 1.1))
  dev.off()
  
  # QC-calculate the gene number of each cell 
  gene_counts_per_cell <- as.numeric(object$nFeature_RNA)
  # apply(object@assays$RNA@counts, 2, function(x) sum(x > 
  #   0))
  # QC-statistics of gene counts per cell
  message("QC-calculate the gene number of each cell")
  print(summary(gene_counts_per_cell))
  gene_counts_summary <- summary(gene_counts_per_cell)
  capture.output(gene_counts_summary, file = paste0(path, "/1.QC/", samplename[i], 
                                                    "/gene_counts_statitics.txt"))
  
  p <- hist(gene_counts_per_cell, breaks = 200)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of gene counts.pdf"))
  # QC-histogram of gene counts distribution
  hist(gene_counts_per_cell, breaks = 200, main = "Distribution of gene counts", 
       xlab = "Gene counts", ylab = "Cell number", col = "blue", xlim = c(0, max(p$breaks) * 
                                                                            1.1), ylim = c(0, max(p$counts) * 1.1))
  dev.off()
  rm(gene_counts_per_cell)
  # QC-visualize UMI counts, gene counts, mitochondria pct, HB pct individually
  gc()
  p <- VlnPlot(object, col = "blue", features = c("nFeature_RNA", "nCount_RNA", 
                                                  "percent.mt", "percent.HB"), ncol = 4, pt.size = 0.01)
  print(p)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct, HB pct.pdf"))
  print(p)
  dev.off()
  
  
  plot1 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot2 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.mt", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot3 <- FeatureScatter(object, feature1 = "nCount_RNA", feature2 = "percent.HB", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot4 <- FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.mt", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  plot5 <- FeatureScatter(object, feature1 = "nFeature_RNA", feature2 = "percent.HB", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct, HB pct in feature scatter.pdf"))
  p <- CombinePlots(plots = list(plot1, plot2, plot3, plot4, plot5), ncol = 3, 
                    legend = "none")
  print(p)
  dev.off()
  print(p)
  nmads=3
  
  
  set.seed(100)
  # Setting up the parameters for consistency with denoisePCA();
  # this can be changed depending on your feature selection scheme.
  dbl.dens <- computeDoubletDensity(object@assays$RNA@counts)
  sce <- as.SingleCellExperiment(object)
  sce <- scDblFinder(sce)
  table(sce$scDblFinder.class)
  # summary(dbl.dens[which(dbl.dens<quantile(dbl.dens,0.95))])
  # table(dbl.dens>quantile(dbl.dens,0.95))
  object$DoubletScore <- sce$scDblFinder.score
  object$scDblFinder.class <- ifelse(sce$scDblFinder.class == 'doublet',T,F)
  rm(sce)
  # table(dbl.dens>quantile(dbl.dens,0.90),object@active.ident)
  p <- hist(dbl.dens, breaks = 200)
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of DoubletScore.pdf"))
  # QC-histogram of gene counts distribution
  hist(dbl.dens, breaks = 200, main = "Distribution of DoubletScore", 
       xlab = "Gene counts", ylab = "Cell number", col = "blue", xlim = c(0, max(p$breaks) * 
                                                                            1.1), ylim = c(0, max(p$counts) * 1.1))
  dev.off()
  
  pc=ifelse(ncol(object)>10000,0.924,
            ifelse(ncol(object)>8000,0.939,
                   ifelse(ncol(object)>6000,0.954,
                          ifelse(ncol(object)>4000,0.969,
                                 ifelse(ncol(object)>2000,0.984,                                        ifelse(ncol(object)>1000,0.992,0.996))))))
  
  
  libsize.drop <- object$nCount_RNA < 500 | object$nCount_RNA > quantile(object$nCount_RNA,pc)
  feature.drop <- object$nFeature_RNA > quantile(object$nFeature_RNA,pc)
  mt.drop <- object$percent.mt > 20
  DoubletS <- object$scDblFinder.class
  keep <- !(libsize.drop | feature.drop | mt.drop | DoubletS)
  filter <- data.frame(ByLibSize = sum(libsize.drop), ByFeature = sum(feature.drop), ByMT = sum(mt.drop), ByDoubletScore=sum(DoubletS),
                       Remaining = sum(keep))
  write.csv(filter,file = paste0(path,"/1.QC/", samplename[i], "/cell filtration.csv"))
  
  # 
  message("Subset seurat object: ", samplename[i])
  MTloc <- c(MTloc,grep('^RP[SL]',rownames(object@assays$RNA@counts)))
  datanew <- object@assays$RNA@counts[-MTloc, keep]
  objectnew <- CreateSeuratObject(counts = datanew, project = samplename[i], min.cells = 3, min.features = 200,)
  object.list[[samplename[i]]] <- objectnew
  
  # object.list[[samplename[i]]] <- object[-MTloc, keep]
  
  message("dim of subset object:", dim(object.list[[samplename[i]]])[1],' ', dim(object.list[[samplename[i]]])[2])
  
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts, mitochondria pct, HB pct for newseurat.pdf"))
  p <- VlnPlot(object.list[[samplename[i]]], col = "blue", features = c("nFeature_RNA", 
                                                                        "nCount_RNA"), ncol = 2, pt.size = 0.01)
  print(p)
  dev.off()
  print(p)
  plot1 <- FeatureScatter(object.list[[samplename[i]]], feature1 = "nCount_RNA", feature2 = "nFeature_RNA", 
                          cols = "black", pt.size = 0.01) + NoLegend()
  
  pdf(file = paste0(path, "/1.QC/", samplename[i], "/Distribution of UMI counts, gene counts in feature scatter after.pdf"))
  p <- plot1
  print(p)
  dev.off()
  print(p)
  rm(object)
  gc()
  message("creat seurat object ", samplename[i], " finished!")
  Sys.sleep(1/1000)
  pb$tick()
}

save(object.list,file = 'object.list.rda')

dir.create('singlecell')
dir.create('singlecell/0.rawdata')

dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
dirlist <- dir(path = 'singlecell/0.rawdata',full.names = T)
path='singlecell'
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dirlist[which(samplename %in% clinical$X.Sample_geo_accession)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]

object.list <- lapply(X = object.list, FUN = function(object) { 
  samplename = as.character(object@meta.data$orig.ident[1])
  # print(samplename)
  if(grepl('-1',colnames(object@assays$RNA@counts)[1])){
    object <- RenameCells(object,new.names = str_replace_all(colnames(object@assays$RNA@counts), '-1',paste0('-',samplename)))
  }else{
    object <- RenameCells(object,new.names = paste0(colnames(object@assays$RNA@counts), '-',samplename))       
  }
  object <- NormalizeData(object) 
  object <- FindVariableFeatures(object, selection.method = "vst", nfeatures = 2000) 
  features <- object@assays$RNA@var.features
  write.csv(features,file=paste0('singlecell','/3.integrate_and_cluster/1.visualize_high_varible_genes/',samplename,'top2000.csv'))
  return(object)
})

features <- SelectIntegrationFeatures(object.list = object.list)
ifnb.anchors <- FindIntegrationAnchors(object.list = object.list, anchor.features = features)
save(ifnb.anchors,file=paste0('Step3.ifnb.anchors.rda'))

pb <- progress_bar$new(
  format = "  over [:bar] :percent time :elapsed",
  total = length(object.list), clear = FALSE, width= 60)
for(i in 1:length(object.list)){
  object.subset <- object.list[[i]]
  samplename = as.character(object.subset@meta.data$orig.ident[1])
  message("visualize high varible genes: ",samplename)
  if (!file.exists(paste0('singlecell','/3.integrate_and_cluster/','1.visualize_high_varible_genes'))){
    dir.create(paste0('singlecell','/3.integrate_and_cluster/','1.visualize_high_varible_genes'))
  }
  top20 <- head(VariableFeatures(object.subset), 20)
  plot6 <- VariableFeaturePlot(object.subset, pt.size = 2, cols = c("black", "blue"))
  LabelPoints(plot = plot6, points = top20, repel = TRUE)
  plot7 <- LabelPoints(plot = plot6, points = top20, repel = TRUE)+ NoLegend()
  pdf(paste0('singlecell','/3.integrate_and_cluster/','1.visualize_high_varible_genes/',samplename," visualize high varible genes.pdf"),width=17,height =9)
  p <- CombinePlots(plots = list(plot6, plot7), legend = "bottom")
  print(p)
  dev.off()
  print(p)
  Sys.sleep(1 / 1000)
  pb$tick()
}

load('Step3.ifnb.anchors.rda')
PBMC.combined <- IntegrateData(anchorset = ifnb.anchors)
DefaultAssay(PBMC.combined) <- "integrated"

load('Step3.PBMC.combined.rda')
dir.create('singlecell/4.normalization')
dir.create('singlecell/4.normalization/2.PCA_analysis')

PBMC.combined <- ScaleData(PBMC.combined, verbose = FALSE,vars.to.regress = c("nCount_RNA",'nFeature_RNA'))

PBMC.combined <- RunPCA(PBMC.combined, npcs = 50, verbose = FALSE, do.print = TRUE, 
                        pcs.print = 1:5, 
                        genes.print = 5)
object <- PBMC.combined
save(PBMC.combined,file='Step6.PBMC.combined.rda')

pdf(paste0('singlecell', "/4.normalization/", "2.PCA_analysis/10.PCA gene（top9）.pdf"), width = 15, height = 20)
print(VizDimLoadings(object, dims = 1:9, reduction = "pca",ncol = 3))
dev.off()


pdf(file = paste0('singlecell', "/4.normalization/", "2.PCA_analysis", "/11.PCA（top9）.pdf"), 
    width = 9, height = 6)
DimHeatmap(object, dims = 1:9, cells = 100, reduction = "pca", ncol = 3, balanced = TRUE)
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "2.PCA_analysis", "/12.PCA.pdf"), 
    width = 8, height = 6)
print(DimPlot(object, dims = 1:2, reduction = "pca",group.by = 'orig.ident') )
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "2.PCA_analysis", "/13.Examine and visualize PCA singlecells with ElbowPlot.pdf"), 
    width = 9, height = 9)
print(ElbowPlot(object, ndims = 50, reduction = "pca"))
dev.off()

rm(list = ls())
load(file=paste0('Step6.PBMC.combined.rda'))
# object <- PBMC.combined
save.fig=F
# 
dir.create(paste0('singlecell', "/4.normalization/3.cluster"))
PBMC.combined$orig.ident <- factor(PBMC.combined$orig.ident)

PBMC.combined <- FindNeighbors(PBMC.combined, reduction = "pca", dims = 1:30)

PBMC.combined <- FindClusters(
  object = PBMC.combined,
  reduction = "pca",
  resolution = 1,
  dims= 1:30,
  save.SNN = TRUE
)
save(PBMC.combined,file=paste0('Step5.FindNeighbors.resolution.PBMC.combined.rda'))



resolution=1
integrated_snn_res <- paste0("integrated_snn_res.", resolution)
Idents(PBMC.combined) <- PBMC.combined@meta.data[[integrated_snn_res]]
#run non linear dimensionality reduction-UMAP
PBMC.combined <- RunUMAP(PBMC.combined, reduction = "pca", dims = 1:30)
DimPlot(PBMC.combined, reduction = "umap", label = TRUE)

PBMC.combined <- RunTSNE(PBMC.combined, dims = 1:30)
DimPlot(PBMC.combined, reduction = "tsne", label = TRUE)
save(PBMC.combined,file = 'Step6.PBMC.combined.rda')
object <- PBMC.combined
pdf(file = paste0('singlecell', "/4.normalization/", "3.cluster", "/14. tsne.pdf"), width = 8, height = 6)
print(DimPlot(object, reduction = "tsne", label = TRUE))
dev.off()
p1 <- DimPlot(object, reduction = "tsne",  pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, label = TRUE)
pdf(file = paste0('singlecell', "/4.normalization/", "3.cluster", "/15. tsne.pdf"), width = 8, height = 6)
print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
dev.off()

pdf(file = paste0('singlecell', "/4.normalization/", "3.cluster", "/16. umap.pdf"), width = 8, height = 6)
print(DimPlot(object, reduction = "umap", label = TRUE))
dev.off()
# Visualization

# Patients <- ifelse(grepl('HCCP',colnames(object)),'HCC Patients','Normal Liver Patients')
# names(Patients) <- colnames(object)
# object[['Patients']] <- Patients
p1 <- DimPlot(object, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "umap", label = TRUE)
pdf(file = paste0('singlecell', "/4.normalization/", "3.cluster", "/17. umap.pdf"), width = 8, height = 6)
print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
dev.off()

object.all.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.6, logfc.threshold = 0.25,test.use = 'roc')
object.all.markers <- object.all.markers[order(object.all.markers$cluster,object.all.markers$pct_diff,decreasing = T),]

object.all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = pct_diff) %>%
  kable("html",caption = '<center>**10. marker genes top9**</center>') %>%
  kable_styling(bootstrap_options = "striped", full_width = F)%>%
  scroll_box(width = "1000px", height = "500px")

write.csv(object.all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = pct_diff),file=paste0('singlecell','/4.normalization/4.cluster_markers/10. marker genes top9.csv'))

identtotal = sort(as.numeric(as.character(unique(object@active.ident))))
for (ident in identtotal) {
  # message(paste0("findmarkers of ", ident))
  # findmarkers of each cluster
  # cluster.markers <- FindMarkers(object, ident.1 = ident, only.pos = TRUE,
  #     min.pct = 0.25, logfc.threshold = 0.7)
  cluster.markers <- object.all.markers[which(object.all.markers$cluster == ident),]
  # print(head(cluster.markers, n = 10))
  
  # cell marker gene
  # matchcelltypegenebelongto <- ConfirmCelltypebyCellmarker(cluster.markers, cell_marker = cell_marker)
  features = as.character(cluster.markers$gene[1:9])
  # visualization of cell marker-violin plot
  pdf(file = paste0('singlecell', "/4.normalization/", "4.cluster_markers", "/20. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(VlnPlot(object, features = features, log = TRUE, ncol = 3))
  dev.off()
  
  # visualization of cell marker-ridge plot
  pdf(file = paste0('singlecell', "/4.normalization/", "4.cluster_markers", "/20. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(RidgePlot(object, features = features, ncol = 3))
  dev.off()
  
  # visualization of cell marker-Dot plot
  pdf(file = paste0('singlecell', "/4.normalization/", "4.cluster_markers", "/20. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(DotPlot(object, features = features, cols = c("lightgrey", 
                                                      "blue")) + RotatedAxis())
  dev.off()
  # visualization of cell marker-fetureplot/umap
  # print(FeaturePlot(object = object, reduction = "umap", features = as.character(cluster.markers$gene[1:9]), 
  # pt.size = 0.01))
  
  # visualization of cell marker-fetureplot/tsne
  pdf(file = paste0('singlecell', "/4.normalization/", "4.cluster_markers", "/20. cluster ",ident," top9 marker tsne.pdf"), width = 10, height = 10)
  print(FeaturePlot(object = object, reduction = "tsne", features = features, 
                    pt.size = 0.01, ncol = 3))
  dev.off()
}
ident=1
cluster.markers <- object.all.markers[which(object.all.markers$cluster == ident),]
# print(head(cluster.markers, n = 10))

# cell marker gene
# matchcelltypegenebelongto <- ConfirmCelltypebyCellmarker(cluster.markers, cell_marker = cell_marker)

# visualization of cell marker-violin plot
VlnPlot(object, features = as.character(cluster.markers$gene[1:9]), log = TRUE, 
        ncol = 3)

# visualization of cell marker-ridge plot
RidgePlot(object, features = as.character(cluster.markers$gene[1:9]), ncol = 3)

# visualization of cell marker-Dot plot
DotPlot(object, features = as.character(cluster.markers$gene[1:9]), cols = c("lightgrey", 
                                                                             "blue")) + RotatedAxis()

# visualization of cell marker-fetureplot/umap
# print(FeaturePlot(object = object, reduction = "umap", features = as.character(cluster.markers$gene[1:9]), 
# pt.size = 0.01))

# visualization of cell marker-fetureplot/tsne
FeaturePlot(object = object, reduction = "tsne", features = as.character(cluster.markers$gene[1:9]), 
            pt.size = 0.01, ncol = 3)
rm(list = ls())
load(file=paste0('Step6.PBMC.combined.rda'))
save.fig=F
#
dir.create(paste0('singlecell', "/4.normalization/",'5.cluster_named'))
symbol.trans <- read.table('symbol_trans.txt',header = F,stringsAsFactors = F)
symbol.trans.fun <- function(features){
  loc <- match(features,symbol.trans[,1])
  match <- symbol.trans[loc[which(!is.na(loc))],2]
  features[which(!is.na(loc))] <- match
  return(features)
}

load('Step6.PBMC.combined.rda')
resolution=1
integrated_snn_res <- paste0("integrated_snn_res.", resolution)
Idents(PBMC.combined) <- PBMC.combined@meta.data[[integrated_snn_res]]

object <- PBMC.combined
# features <- c('PTPRC','CD7','EPCAM','PDGFRA','PDGFRB')
features <- unique(c(
  'PECAM1', 'PLVAP', 'VWF', 'CDH5', #ENCs
  'COL1A2','COL1A1','DCN', 'LUM', #Fibroblast cells
  'EPCAM','KRT7','KRT8','KRT19','CDH1', #EPCs
  'CD79A','MS4A1','VPREB3','CD37', #B Cells
  'FCER1G','C1QA','C1QB', 'AIF1', #Myeloid cells
  'GNLY','NKG7','KLRD1','GZMB', #NK Cells
  'TRBC2','CD2','CD3D','CD3E','CD3G','CD4','CD8A' # T cells
))
features <- symbol.trans.fun(features)
object@meta.data$ident <- object@active.ident
if(save.fig == T){
  pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/42. marker.pdf"), width = 12, height = 6)
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
object <- PBMC.combined
clustername <- c('T cells','T cells','NK cells','T cells','Myeloid cells','B cells','Myeloid cells','EPCs','T cells','T cells','T cells','T cells','T cells','Myeloid cells','Myeloid cells','Fibroblast cells','Myeloid cells','Myeloid cells','T cells','Myeloid cells','ENCs','T cells','Myeloid cells','EPCs','T cells')

new.cluster.ids <- clustername
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)

# features <- c('PTPRC','CD7','EPCAM','PDGFRA','PDGFRB')
features <- unique(c('GNG11','IFI27','TM4SF1','CLDN5','PECAM1', 'PLVAP', 'VWF', 'CDH5','TAGLN','ACTA2','TPM2','MYL9','BGN','FHL2','AEBP1','COL1A2','COL1A1','DCN', 'LUM','EPCAM','KRT7','KRT8','KRT19','CDH1','CD79A','MS4A1','CD22',"IGLC2","IGKC","IGLC3","IGHA1",'VPREB3','CD37','FCER1G','CD14','CXCL8','C1QA','C1QB', 'AIF1','ACTA2','GNLY','CCL3','NKG7','KLRD1','GZMB','GNLY','NKG7','SORCS1', 'MYOT','EPCAM','KRT8','MUC2','IL7R','TRBC2','CD2','CD3D','CD3E','CD3G','CD4','CD8A'))
features <- symbol.trans.fun(features)
object@meta.data$ident <- object@active.ident
if(save.fig == T){
  pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/43. marker.pdf"), width = 12, height = 6)
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
# PBMC.combined <- object
# save(PBMC.combined,file = 'Step6.PBMC.combined.rda')
object <- PBMC.combined
clustername <- c('T cells','T cells','NK cells','T cells','Myeloid cells','B cells','Myeloid cells','EPCs','T cells','T cells','T cells','T cells','T cells','Myeloid cells','Myeloid cells','Fibroblast cells','Myeloid cells','Myeloid cells','T cells','Myeloid cells','ENCs','T cells','Myeloid cells','EPCs','T cells')

new.cluster.ids <- clustername
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)

features <- unique(c(
  'PECAM1', 'PLVAP', 'VWF', 'CDH5', #ENCs
  'COL1A2','COL1A1','DCN', 'LUM', #Fibroblast cells
  'EPCAM','KRT7','KRT8','KRT19','CDH1', #EPCs
  'CD79A','MS4A1','VPREB3','CD37', #B Cells
  'FCER1G','C1QA','C1QB', 'AIF1', #Myeloid cells
  'GNLY','NKG7','KLRD1','GZMB', #NK Cells
  'TRBC2','CD2','CD3D','CD3E','CD3G','CD4','CD8A' # T cells
))
features <- symbol.trans.fun(features)
object@meta.data$ident <- object@active.ident
if(save.fig == T){
  pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/43. marker.pdf"), width = 12, height = 6)
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
# PBMC.combined <- object
# save(PBMC.combined,file = 'Step6.PBMC.combined.rda')
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/44. tsne.pdf"), width = 6, height = 6)
print(DimPlot(object, reduction = "tsne", label = TRUE, pt.size = 0.001) + NoLegend())
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/45.PCA.pdf"), 
    width = 8, height = 6)
print(DimPlot(object, dims = 1:2, reduction = "pca") )
dev.off()
p1 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, label = TRUE)

pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/46. tsne.pdf"), width = 8, height = 5)
print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/47. tsne.pdf"), width = 10, height = 12)
print(DimPlot(object, reduction = "tsne", split.by = "orig.ident", ncol = 2))
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/25. umap.pdf"), width = 6, height = 6)
print(DimPlot(object, reduction = "umap", label = TRUE, pt.size = 0.001) + NoLegend())
dev.off()
p1 <- DimPlot(object, reduction = "umap", pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "umap", pt.size = 0.001, label = TRUE)

pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/26. umap.pdf"), width = 8, height = 5)
print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
dev.off()
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]
samplenamenormal <- samplename[which(group %in%'Normal')]
samplenametumor <- samplename[which(group %in%'Tumor')]
PBMC.combined <- object
a <- table(PBMC.combined$orig.ident,PBMC.combined@active.ident)

a <- round(a/rowSums(a),4)
Cellratio <- as.data.frame(a)
Cellratio$group <- ifelse(Cellratio$Var1 %in% samplenamenormal,'Normal','Tumor')
Cellratio$group <- factor(Cellratio$group,levels = c('Tumor','Normal'))
Cellratio$Var1 <- factor(Cellratio$Var1,levels = c(samplenamenormal,samplenametumor))
library(ggpubr)
p <- ggboxplot(Cellratio, x = "Var2", y = "Freq", 
               color = "group", ylab = "cell ratio", 
               xlab = "", add = "jitter", size = 1, axis.line = 2) + 
  stat_compare_means(label = "p.signif", label.y = max(Cellratio[, 
                                                                 3]), aes(group = group),method = 'kruskal.test') + rotate_x_text(60)

if(save.fig == T){
  pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/27-1. .pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}else{
  print(p)
}
# Cellratio <- Cellratio[order(Cellratio$group),]
# Cellratio$Var1 <- factor(Cellratio$Var1,levels = unique(Cellratio$Var1),ordered = T)
p1 <- ggplot(data = Cellratio,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat="identity",position = "fill")+
  ##
  geom_hline(yintercept=0.25,linetype=2,size=1)+
  geom_hline(yintercept=0.50,linetype=2,size=1)+
  geom_hline(yintercept=0.75,linetype=2,size=1)+
  ##
  # coord_flip()+
  xlab("")+
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  theme(
    ##
    legend.position = "bottom",legend.title = element_blank(),
    legend.text=element_text(size = 14,colour = "black"),
    ####
    line = element_line(colour = "black", size = 1, linetype = 1), 
    ###
    axis.title=element_text(size = 16,face="bold",colour = "black"),
    ###
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.text.x = element_text(size = 14,colour = "black",angle = 45,hjust = 1))

sample=c(samplenamenormal,samplenametumor)
annotation_col<- data.frame(sample ,
                            Var1 =ifelse(sample %in% samplenamenormal,'Normal','Tumor'), 
                            Var2=c(1:length(unique(Cellratio$Var1))),
                            value=1)
annotation_col$Var1 <- factor(annotation_col$Var1,levels = c('Tumor','Normal'))
p2 <- ggplot()+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+scale_fill_brewer(palette = "Set1")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  # scale_fill_manual(labels=c("preoperative","postoperative_1M","postoperative_3M","postoperative_6M"),
  #                   values =c("preoperative"="#F8766D",
  #                             "postoperative_1M"="#7CAE00",
  #                             "postoperative_3M"="#00BFC4",
  #                             "postoperative_6M"="#C77CFF"))+
  theme(plot.margin = ggplot2::margin(1, 1, 0, 0, "cm"))
# p2
library(cowplot)
#
if(save.fig == T){
  pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/27-2. rate.pdf"), width = 12, height = 6)
  print(plot_grid(p2,p1,
                  rel_heights = c(0.5,1),
                  label_x=0,
                  label_y=1,
                  scale = c(1.069, 1),
                  align = 'v',ncol = 1,greedy = F) )
  dev.off()
}else{
  print(plot_grid(p2,p1,
                  rel_heights = c(0.5,1),
                  label_x=0,
                  label_y=1,
                  scale = c(1.069, 1),
                  align = 'v',ncol = 1,greedy = F) )
}
features <- c('BMP8B','DENND3','FAM155B','IL22RA2','ITM2A','MAP3K8','OTX1','SCARF1','TPO','ZMYND15','ABCB4','CBFA2T3','FCAMR','ISL2','MC4R','TAS1R3')
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]
samplenamenormal <- samplename[which(group %in%'Normal')]
samplenametumor <- samplename[which(group %in%'Tumor')]
ifelse(Cellratio$Var1 %in% samplenamenormal,'Normal','Tumor')
PBMC.combined$group <- ifelse(PBMC.combined$orig.ident %in% samplenamenormal,'Normal','Tumor')
PBMC.combined@meta.data$ident <- paste(PBMC.combined@active.ident,PBMC.combined$group,sep = ' | ')
if(save.fig == T){
  pdf(file = paste0( 'singlecell', "/4.normalization/", "5.cluster_named", "/28. 16.pdf"), width = 10, height = 6)
  print(DotPlot(PBMC.combined, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(PBMC.combined, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
data <- PBMC.combined@assays$RNA@data
genes <- intersect(rownames(data), features)
boxplotdata <- t(data[genes, ])
boxdatalist <- lapply(colnames(boxplotdata), function(x) data.frame(expression = log2(boxplotdata[, x] + 1), 
                                                                    group = PBMC.combined@meta.data$group, 
                                                                    gene = Idents(PBMC.combined)))
# boxdata <- do.call(rbind, boxdata)
i=1
boxdata <- boxdatalist[[i]]
boxdata$group <- factor(boxdata$group, levels = c("Tumor",
                                                  "Normal"))
p = ggboxplot(boxdata, x = "gene", y = "expression", 
              color = "group", ylab = "Gene expression", 
              xlab = "", add = "jitter", size = 1, axis.line = 2) + 
  stat_compare_means(label = "p.signif", label.y = max(boxdata[, 
                                                               1]), aes(group = group)) + rotate_x_text(60)
if(save.fig == T){
  pdf(file = paste0( 'singlecell', "/4.normalization/", "5.cluster_named", "/29. ",colnames(boxplotdata)[i],"boxplots.pdf"), width = 6, height = 6)
  print(p)
  dev.off()
}else{
  print(p)
}

for (i in 1:length(colnames(boxplotdata))) {
  boxdata <- boxdatalist[[i]]
  boxdata$group <- factor(boxdata$group, levels = c("Tumor",
                                                    "Normal"))
  p = ggboxplot(boxdata, x = "gene", y = "expression", 
                color = "group", ylab = "Gene expression", 
                xlab = "", add = "jitter", size = 1, axis.line = 2) + 
    stat_compare_means(label = "p.signif", label.y = max(boxdata[, 
                                                                 1]), aes(group = group)) + rotate_x_text(60)
  if(save.fig == T){
    pdf(file = paste0( 'singlecell', "/4.normalization/", "5.cluster_named", "/29. ",colnames(boxplotdata)[i],"boxplots.pdf"), width = 6, height = 6)
    print(p)
    dev.off()
  }else{
    print(p)
  }
}
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/30.violin1-8.pdf"), width = 8, height = 15)
print(VlnPlot(object, features = features[1:8], log = TRUE, ncol = 1,pt.size = 0))
dev.off()
pdf(file = paste0('singlecell', "/4.normalization/", "5.cluster_named", "/30.9-16.pdf"), width = 8, height = 15)
print(VlnPlot(object, features = features[9:16], log = TRUE, ncol = 1,pt.size = 0))
dev.off()
load('Step6.PBMC.combined.rda')
table(Idents(PBMC.combined))
PBMC.combined.Myeloid <- PBMC.combined[,Idents(PBMC.combined)=='Myeloid cells']
rm(PBMC.combined)
PBMC.combined <- PBMC.combined.Myeloid
rm(PBMC.combined.Myeloid)
dim(PBMC.combined)
DefaultAssay(PBMC.combined) <- 'integrated'
unique(Idents(PBMC.combined))
save(PBMC.combined,file = 'Step10.PBMC.combined.Myeloid.rda')

load('Step10.PBMC.combined.Myeloid.rda')
object <- PBMC.combined
object <- FindVariableFeatures(object,  selection.method = "vst",nfeatures = 1000,verbose = FALSE) 
dir.create(paste0('singlecell', "/5.Myeloidsubclustering"))
dir.create(paste0('singlecell', "/5.Myeloidsubclustering/1.visualize_high_varible_genes"))
top20 <- head(VariableFeatures(object), 20)
plot6 <- VariableFeaturePlot(object, pt.size = 2, cols = c("black", "blue"))
# LabelPoints(plot = plot6, points = top20, repel = TRUE)
plot7 <- LabelPoints(plot = plot6, points = top20, repel = TRUE) + NoLegend()
if(save.fig == T){
  pdf(paste0('singlecell', "/5.Myeloidsubclustering/", "1.visualize_high_varible_genes/", 
             "272. visualize high varible genes.pdf"), width = 8, height = 6)
  print(CombinePlots(plots = list(plot6, plot7), legend = "bottom"))
  dev.off()
}else{
  CombinePlots(plots = list(plot6, plot7), legend = "bottom")
}
dir.create('singlecell/5.Myeloidsubclustering/2.PCA_analysis')
object <- RunPCA(object, npcs = 50, verbose = FALSE, do.print = TRUE, 
                 pcs.print = 1:5, genes.print = 5)
pdf(paste0('singlecell', "/5.Myeloidsubclustering/", "2.PCA_analysis/272.PCAgene（top9）.pdf"), width = 15, height = 20)
print(VizDimLoadings(object, dims = 1:9, reduction = "pca",ncol = 3))
dev.off()
# explore the primary sources of heterogeneity in a dataset all PC
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "2.PCA_analysis", "/273.PCA（top9）.pdf"), 
      width = 9, height = 6)
  DimHeatmap(object, dims = 1:9, cells = 100, reduction = "pca", ncol = 3, balanced = TRUE)
  dev.off()
}else{
  DimHeatmap(object, dims = 1:9, cells = 100, reduction = "pca", ncol = 3, balanced = TRUE)
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "2.PCA_analysis", "/274.PCA.pdf"), 
      width = 8, height = 6)
  print(DimPlot(object, dims = 1:2, reduction = "pca",group.by = 'orig.ident') )
  dev.off()
}else{
  DimPlot(object, dims = 1:2, reduction = "pca",group.by = 'orig.ident') 
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "2.PCA_analysis", "/275.Examine and visualize PCA results with ElbowPlot.pdf"), 
      width = 9, height = 9)
  print(ElbowPlot(object, ndims = 50, reduction = "pca"))
  dev.off()
}else{
  ElbowPlot(object, ndims = 50, reduction = "pca")
}
dir.create(paste0('singlecell', "/5.Myeloidsubclustering/3.cluster"))
PBMC.combined$orig.ident <- factor(PBMC.combined$orig.ident)
object <- FindNeighbors(object, reduction = "pca", dims = 1:30)
object <- FindClusters(
  object = object,
  reduction = "pca",
  resolution = 1,
  dims= 1:30,
  save.SNN = TRUE
)
pdf(paste0('singlecell', "/5.Myeloidsubclustering/3.cluster/clustree.pdf"),width = 8,height = 6)
clustree(object@meta.data, prefix = "integrated_snn_res.",prop_filter =0.1)
dev.off()
resolution=1
integrated_snn_res <- paste0("integrated_snn_res.", resolution)
Idents(object) <- object@meta.data[[integrated_snn_res]]
#run non linear dimensionality reduction-UMAP
object <- RunUMAP(object, reduction = "pca", dims = 1:30)
DimPlot(object, reduction = "umap", label = TRUE)

object <- RunTSNE(object, dims = 1:30,perplexity = 40)
DimPlot(object, reduction = "tsne", label = TRUE)
PBMC.combined <- object
save(PBMC.combined,file = 'Step10.PBMC.combined.Myeloid.filter.rda')
pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/276. tsne.pdf"), width = 8, height = 6)
print(DimPlot(object, reduction = "tsne", label = TRUE))
dev.off()
# Visualization
p1 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, label = TRUE)
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/277. tsne.pdf"), width = 8, height = 6)
  print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
  dev.off()
}else{
  CombinePlots(plots = list(p1, p2), legend = "bottom")
}

if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/278. umap.pdf"), width = 8, height = 6)
  print(DimPlot(object, reduction = "umap", label = TRUE))
  dev.off()
}else{
  DimPlot(object, reduction = "umap", label = TRUE)
}
# Visualization

# Patients <- ifelse(grepl('HCCP',colnames(object)),'HCC Patients','Normal Liver Patients')
# names(Patients) <- colnames(object)
# object[['Patients']] <- Patients
p1 <- DimPlot(object, reduction = "umap", group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "umap", label = TRUE)
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/279. umap.pdf"), width = 8, height = 6)
  print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
  dev.off()
}else{
  CombinePlots(plots = list(p1, p2), legend = "bottom")
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/280. tsne.pdf"), width = 10, height = 12)
  print(DimPlot(object, reduction = "tsne", split.by = "orig.ident", ncol = 4))
  dev.off()
}else{
  DimPlot(object, reduction = "tsne", split.by = "orig.ident", ncol = 4)
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "3.cluster", "/281. umap.pdf"), width = 10, height = 12)
  print(DimPlot(object, reduction = "umap", split.by = "orig.ident", ncol = 4))
  dev.off()
}else{
  DimPlot(object, reduction = "umap", split.by = "orig.ident", ncol = 3)
}
dir.create(paste0('singlecell', "/5.Myeloidsubclustering/",'4.cluster_markers'))
PBMC.combined.all.markers <- FindAllMarkers(object, only.pos = TRUE, min.pct = 0.5, logfc.threshold = 0.25,test.use = 'roc')
# PBMC.combined.all.markers %>% group_by(cluster) %>% top_n(n = 2, wt = avg_diff)

#write the marker genes to the txt file
write.table(PBMC.combined.all.markers,file=paste0('singlecell', "/5.Myeloidsubclustering/",'4.cluster_markers',"/markers_allcluster.txt"),row.names = T,quote = FALSE,sep = '\t')

object.all.markers <- read.table(paste0('singlecell', "/5.Myeloidsubclustering/", 
                                        "4.cluster_markers", "/markers_allcluster.txt"),header = T)
object.all.markers$pct_diff <- object.all.markers$pct.1 - object.all.markers$pct.2
object.all.markers <- object.all.markers[order(object.all.markers$cluster,object.all.markers$pct_diff,decreasing = T),]

object.all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = pct_diff) %>%
  kable("html",caption = '<center>**110. marker genes top9**</center>') %>%
  kable_styling(bootstrap_options = "striped", full_width = F)%>%
  scroll_box(width = "1000px", height = "500px")

write.csv(object.all.markers %>% group_by(cluster) %>% top_n(n = 9, wt = pct_diff),file=paste0('singlecell','/5.Myeloidsubclustering/4.cluster_markers/110. marker genes top9.csv'))

# findmarkers of each cluster
identtotal = sort(as.numeric(as.character(unique(object@active.ident))))
for (ident in identtotal) {
  # message(paste0("findmarkers of ", ident))
  # findmarkers of each cluster
  # cluster.markers <- FindMarkers(object, ident.1 = ident, only.pos = TRUE,
  #     min.pct = 0.25, logfc.threshold = 0.7)
  cluster.markers <- object.all.markers[which(object.all.markers$cluster == ident),]
  # print(head(cluster.markers, n = 10))
  if(nrow(cluster.markers) == 0) next
  # 查看cell marker gene是那些组织的标志物
  # matchcelltypegenebelongto <- ConfirmCelltypebyCellmarker(cluster.markers, cell_marker = cell_marker)
  features = as.character(cluster.markers$gene[1:9])
  # visualization of cell marker-violin plot
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "4.cluster_markers", "/282. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(VlnPlot(object, features = features, log = TRUE, ncol = 3))
  dev.off()
  
  # visualization of cell marker-ridge plot
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "4.cluster_markers", "/282. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(RidgePlot(object, features = features, ncol = 3))
  dev.off()
  
  # visualization of cell marker-Dot plot
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "4.cluster_markers", "/282. cluster ",ident," top9 marker.pdf"), width = 10, height = 10)
  print(DotPlot(object, features = features, cols = c("lightgrey", 
                                                      "blue")) + RotatedAxis())
  dev.off()
  # visualization of cell marker-fetureplot/umap
  # print(FeaturePlot(object = object, reduction = "umap", features = as.character(cluster.markers$gene[1:9]), 
  # pt.size = 0.01))
  
  # visualization of cell marker-fetureplot/tsne
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "4.cluster_markers", "/282. cluster ",ident," top9 marker tsne.pdf"), width = 10, height = 10)
  print(FeaturePlot(object = object, reduction = "tsne", features = features, 
                    pt.size = 0.01, ncol = 3))
  dev.off()
}
ident=4
cluster.markers <- object.all.markers[which(object.all.markers$cluster == ident),]
# print(head(cluster.markers, n = 10))

# cell marker gene
# matchcelltypegenebelongto <- ConfirmCelltypebyCellmarker(cluster.markers, cell_marker = cell_marker)

# visualization of cell marker-violin plot
VlnPlot(object, features = as.character(cluster.markers$gene[1:9]), log = TRUE, 
        ncol = 3)

# visualization of cell marker-ridge plot
RidgePlot(object, features = as.character(cluster.markers$gene[1:9]), ncol = 3)

# visualization of cell marker-Dot plot
DotPlot(object, features = as.character(cluster.markers$gene[1:9]), cols = c("lightgrey", 
                                                                             "blue")) + RotatedAxis()

# visualization of cell marker-fetureplot/umap
# print(FeaturePlot(object = object, reduction = "umap", features = as.character(cluster.markers$gene[1:9]), 
# pt.size = 0.01))

# visualization of cell marker-fetureplot/tsne
FeaturePlot(object = object, reduction = "tsne", features = as.character(cluster.markers$gene[1:9]), 
            pt.size = 0.01, ncol = 3)
load('Step10.PBMC.combined.Myeloid.filter.rda')
resolution=1
integrated_snn_res <- paste0("integrated_snn_res.", resolution)
Idents(object) <- object@meta.data[[integrated_snn_res]]
PBMC.combined <- PBMC.combined[,-which(PBMC.combined@active.ident %in% c(11,12))]
save(PBMC.combined,file = 'Step10.PBMC.combined.Myeloid.filter.rda')
save.fig=F
load('Step10.PBMC.combined.Myeloid.filter.rda')
object <- PBMC.combined
resolution=1
integrated_snn_res <- paste0("integrated_snn_res.", resolution)
Idents(object) <- object@meta.data[[integrated_snn_res]]
identtotal = sort(as.numeric(as.character(unique(object@active.ident))))
dir.create('singlecell/5.Myeloidsubclustering/5.cluster_named')
symbol.trans <- read.table('symbol_trans.txt',header = F,stringsAsFactors = F)
symbol.trans.fun <- function(features){
  loc <- match(features,symbol.trans[,1])
  match <- symbol.trans[loc[which(!is.na(loc))],2]
  features[which(!is.na(loc))] <- match
  return(features)
}
object <- PBMC.combined
# features <- c('PTPRC','CD7','EPCAM','PDGFRA','PDGFRB')
features <- unique(c('APOE', 'C1QA','C1QB','S100A8','CXCL8', "G0S2", "IL1B", "EREG", "TIMP1","TPSB2",'TPSAB1','CPA3'))
features <- symbol.trans.fun(features)
object@meta.data$ident <- object@active.ident
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/42. marker.pdf"), width = 12, height = 6)
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
object <- PBMC.combined
clustername <- c('M2-like macrophages','M1-like macrophages','M2-like macrophages','Mast cells','M1-like macrophages','M2-like macrophages','M2-like macrophages','M2-like macrophages','M1-like macrophages','M2-like macrophages','M1-like macrophages')

new.cluster.ids <- clustername
names(new.cluster.ids) <- levels(object)
object <- RenameIdents(object, new.cluster.ids)

# features <- c('PTPRC','CD7','EPCAM','PDGFRA','PDGFRB')
features <- unique(c('APOE', 'C1QA','C1QB','S100A8','CXCL8', "G0S2", "IL1B", "EREG", "TIMP1","TPSB2",'TPSAB1','CPA3'))
features <- symbol.trans.fun(features)
object@meta.data$ident <- object@active.ident
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/43. 类marker.pdf"), width = 12, height = 6)
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(object, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
# PBMC.combined <- object
# save(PBMC.combined,file = 'Step10.PBMC.combined.Myeloid.filter.rda')
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/44.tsne.pdf"), width = 6, height = 6)
  print(DimPlot(object, reduction = "tsne", label = TRUE, pt.size = 0.001) + NoLegend())
  dev.off()
}else{
  DimPlot(object, reduction = "tsne", label = TRUE, pt.size = 0.001) + NoLegend()
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/45.PCA.pdf"), 
      width = 8, height = 6)
  print(DimPlot(object, dims = 1:2, reduction = "pca") )
  dev.off()
}else{
  DimPlot(object, dims = 1:2, reduction = "pca") 
}
p1 <- DimPlot(object, reduction = "tsne",  pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "tsne", pt.size = 0.001, label = TRUE)
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/46. tsne.pdf"), width = 8, height = 5)
  print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
  dev.off()
}else{
  CombinePlots(plots = list(p1, p2), legend = "bottom")
}
pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/47. tsne.pdf"), width = 10, height = 12)
print(DimPlot(object, reduction = "tsne", split.by = "orig.ident", ncol = 2))
dev.off()
# clustername <- c('Monocytes/Macrophage','Mes','Mes','B cells','Monocytes/Macrophage')
# new.cluster.ids <- clustername
# names(new.cluster.ids) <- levels(object)
# object <- RenameIdents(object, new.cluster.ids)
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/25. umap.pdf"), width = 6, height = 6)
  print(DimPlot(object, reduction = "umap", label = TRUE, pt.size = 0.001) + NoLegend())
  dev.off()
}else{
  DimPlot(object, reduction = "umap", label = TRUE, pt.size = 0.001) + NoLegend()
}
p1 <- DimPlot(object, reduction = "umap",  pt.size = 0.001, group.by = "orig.ident")
p2 <- DimPlot(object, reduction = "umap", pt.size = 0.001, label = TRUE)
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/26. umap.pdf"), width = 8, height = 5)
  print(CombinePlots(plots = list(p1, p2), legend = "bottom"))
  dev.off()
}else{
  CombinePlots(plots = list(p1, p2), legend = "bottom")
}
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]
samplenamenormal <- samplename[which(group %in%'Normal')]
samplenametumor <- samplename[which(group %in%'Tumor')]
PBMC.combined <- object
a <- table(PBMC.combined$orig.ident,PBMC.combined@active.ident)

a <- round(a/rowSums(a),4)
Cellratio <- as.data.frame(a)
Cellratio$group <- ifelse(Cellratio$Var1 %in% samplenamenormal,'Normal','Tumor')
Cellratio$group <- factor(Cellratio$group,levels = c('Tumor','Normal'))
Cellratio$Var1 <- factor(Cellratio$Var1,levels = c(samplenamenormal,samplenametumor))
library(ggpubr)
p <- ggboxplot(Cellratio, x = "Var2", y = "Freq", 
               color = "group", ylab = "cell ratio", 
               xlab = "", add = "jitter", size = 1, axis.line = 2) + 
  stat_compare_means(label = "p.signif", label.y = max(Cellratio[, 
                                                                 3]), aes(group = group),method = 'kruskal.test') + rotate_x_text(60)

if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/27-1. rate.pdf"), width = 8, height = 6)
  print(p)
  dev.off()
}else{
  print(p)
}
# Cellratio$group <- factor(Cellratio$group,levels = unique(Cellratio$group),ordered = T)
# Cellratio <- Cellratio[order(Cellratio$group),]
# Cellratio$Var1 <- factor(Cellratio$Var1,levels = unique(Cellratio$Var1),ordered = T)
p1 <- ggplot(data = Cellratio,aes(x=Var1,y=Freq,fill=Var2))+
  geom_bar(stat="identity",position = "fill")+
  ##
  geom_hline(yintercept=0.25,linetype=2,size=1)+
  geom_hline(yintercept=0.50,linetype=2,size=1)+
  geom_hline(yintercept=0.75,linetype=2,size=1)+
  ##
  # coord_flip()+
  xlab("")+
  scale_fill_brewer(palette = "Set3")+
  theme_classic()+
  theme(
    ##
    legend.position = "bottom",legend.title = element_blank(),
    legend.text=element_text(size = 14,colour = "black"),
    ####
    line = element_line(colour = "black", size = 1, linetype = 1), 
    ###
    axis.title=element_text(size = 16,face="bold",colour = "black"),
    ###
    axis.text.y = element_text(size = 14,colour = "black"),
    axis.text.x = element_text(size = 14,colour = "black",angle = 45,hjust = 1))

sample=c(samplenamenormal,samplenametumor)
annotation_col<- data.frame(sample ,
                            Var1 =ifelse(sample %in% samplenamenormal,'Normal','Tumor'), 
                            Var2=c(1:length(unique(Cellratio$Var1))),
                            value=1)
annotation_col$Var1 <- factor(annotation_col$Var1,levels = c('Tumor','Normal'))
p2 <- ggplot()+
  geom_raster(aes(Var2,value,fill = Var1), data= annotation_col)+ 
  labs(x = "", y = "")+scale_fill_brewer(palette = "Set1")+
  theme_void()+
  theme(legend.position = "top", legend.title = element_blank())+
  scale_x_discrete(name="")+
  # scale_fill_manual(labels=c("preoperative","postoperative_1M","postoperative_3M","postoperative_6M"),
  #                   values =c("preoperative"="#F8766D",
  #                             "postoperative_1M"="#7CAE00",
  #                             "postoperative_3M"="#00BFC4",
  #                             "postoperative_6M"="#C77CFF"))+
  theme(plot.margin = ggplot2::margin(1, 1, 0, 0, "cm"))
# p2
library(cowplot)
#
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/27-2. rate.pdf"), width = 8, height = 6)
  print(plot_grid(p2,p1,
                  rel_heights = c(0.5,1),
                  label_x=0,
                  label_y=1,
                  scale = c(1.069, 1),
                  align = 'v',ncol = 1,greedy = F) )
  dev.off()
}else{
  print(plot_grid(p2,p1,
                  rel_heights = c(0.5,1),
                  label_x=0,
                  label_y=1,
                  scale = c(1.069, 1),
                  align = 'v',ncol = 1,greedy = F) )
}
features <- c('DENND3','CBFA2T3')
clinical <- read.csv('GSE123902_clinical.csv')
group <- ifelse(grepl('NORMAL',clinical$X.Sample_title),'Normal','Tumor')
dirlist <- dir(path = 'singlecell/0.rawdata')
samplename <- unlist(stringr::str_split(dirlist,'_MSK'))[seq(1,length(dirlist)*2,2)]
samplename <- samplename[which(samplename %in% clinical$X.Sample_geo_accession)]
samplenamenormal <- samplename[which(group %in%'Normal')]
samplenametumor <- samplename[which(group %in%'Tumor')]
ifelse(Cellratio$Var1 %in% samplenamenormal,'Normal','Tumor')
PBMC.combined$group <- ifelse(PBMC.combined$orig.ident %in% samplenamenormal,'Normal','Tumor')
PBMC.combined@meta.data$ident <- paste(PBMC.combined@active.ident,PBMC.combined$group,sep = ' | ')
if(save.fig == T){
  pdf(file = paste0( 'singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/28. 2gene.pdf"), width = 6, height = 6)
  print(DotPlot(PBMC.combined, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
  dev.off()
}else{
  print(DotPlot(PBMC.combined, group.by = 'ident', features = features,assay='RNA') + RotatedAxis())
}
data <- PBMC.combined@assays$RNA@data
genes <- intersect(rownames(data), features)
boxplotdata <- t(data[genes, ])
boxdatalist <- lapply(colnames(boxplotdata), function(x) data.frame(expression = log2(boxplotdata[, x] + 1), 
                                                                    group = PBMC.combined@meta.data$group, 
                                                                    gene = Idents(PBMC.combined)))
# boxdata <- do.call(rbind, boxdata)
boxdata <- boxdatalist[[1]]
boxdata$group <- factor(boxdata$group, levels = c("Tumor",
                                                  "Normal"))
p = ggboxplot(boxdata, x = "gene", y = "expression", 
              color = "group", ylab = "Gene expression", 
              xlab = "", add = "jitter", size = 1, axis.line = 2) + 
  stat_compare_means(label = "p.signif", label.y = max(boxdata[, 
                                                               1]), aes(group = group)) + rotate_x_text(60)
if(save.fig == T){
  pdf(file = paste0( 'singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/29. ",colnames(boxplotdata)[1],"box.pdf"))
  print(p)
  dev.off()
}else{
  print(p)
}

boxdata <- boxdatalist[[2]]
boxdata$group <- factor(boxdata$group, levels = c("Tumor",
                                                  "Normal"))
p = ggboxplot(boxdata, x = "gene", y = "expression", 
              color = "group", ylab = "Gene expression", 
              xlab = "", add = "jitter", size = 1, axis.line = 2) + 
  stat_compare_means(label = "p.signif", label.y = max(boxdata[, 
                                                               1]), aes(group = group)) + rotate_x_text(60)
if(save.fig == T){
  pdf(file = paste0( 'singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/29. ",colnames(boxplotdata)[2],"box.pdf"))
  print(p)
  dev.off()
}else{
  print(p)
}
DotPlot(PBMC.combined, group.by = 'ident', features = features,assay='RNA') + RotatedAxis()
# features <- makerlist
# features <- c('PHLDA1RASSF6', 'IF127', 'BAMB1', 'TNFRSF19', 'SMOC2', 'NOTUM', 'PROX1')
# features <- c('CMAH','CD90','CD93','CDCP1','CXCR4','EPCR','Erythropoietin receptor','ESAM','EVI1','FLK1','FLK2','FLT3','GATA2','GFI1','MCL1','MYB','PLZF')
features <- unique(symbol.trans.fun(features))
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/29. 2 tsne.pdf"), width = 6, height = 5.5)
  print(FeaturePlot(object = PBMC.combined, reduction = "tsne", features = features, pt.size = 0.01, cols=cols <- c(rgb(211, 211, 211, 50, maxColorValue=255),rgb(252, 187, 161, 50, maxColorValue=255),rgb(251,106,74,50, maxColorValue=255)),min.cutoff='q10',max.cutoff='q90',
                    ncol = 3,split.by='group',order = T))
  dev.off()
}else{
  FeaturePlot(object = PBMC.combined, reduction = "tsne", features = features, pt.size = 1.5, cols=cols <- c(rgb(211, 211, 211, 50, maxColorValue=255),rgb(252, 187, 161, 50, maxColorValue=255),rgb(251,106,74,50, maxColorValue=255)),min.cutoff='q10',max.cutoff='q90',
              ncol = 3,split.by='group')
}
if(save.fig == T){
  pdf(file = paste0('singlecell', "/5.Myeloidsubclustering/", "5.cluster_named", "/30.vio.pdf"), width = 8, height = 8)
  print(VlnPlot(object, features = features[1:2], log = TRUE, ncol = 1,pt.size = 0))
  dev.off()
}else{
  print(VlnPlot(object, features = features[1:2], log = TRUE, ncol = 1,pt.size = 0))
}













