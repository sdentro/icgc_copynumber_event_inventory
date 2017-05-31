args = commandArgs(T)

samplename = args[1]
summary_table = args[2]
consensus_file = args[3]
outdir = args[4]
wgd_file = args[5]

CENTROMERE_FILE = "~/repo/icgc_paper_figures/ucsc_centromere_mid_hg19.txt"

# bb_subclones_file = "battenberg/0040b1b6-b07a-4b6e-90ef-133523eaf412_subclones.txt"
# samplename = "0040b1b6-b07a-4b6e-90ef-133523eaf412"
# summary_table = "~/Documents/Projects/icgc/copynumber_collated_plot/summary_table_combined_annotations.txt"
# consensus_file = "battenberg/0040b1b6-b07a-4b6e-90ef-133523eaf412.consensus.20170119.somatic.cna.annotated.txt"
# outdir = ""
# wgd_file = "~/Documents/Projects/icgc/consensus_subclonal_copynumber/final_run_testing/final_purity_ploidy_patch/combined_consensus_purities_w_sd_allAnno_finalSlimmed_withWGDuncertainty.txt"

summary_table = readr::read_tsv(summary_table)
wgd_anno = readr::read_tsv(wgd_file)

# set this to 2 for not WGD and 4 to WGD samples
is_male = summary_table$inferred_sex[summary_table$samplename==samplename]=="male"
is_wgd = wgd_anno$wgd_status[wgd_anno$samplename==samplename]=="wgd"

if (is_wgd) {
  ploidy = 4
} else {
  ploidy = 2
}

# get centromere mid: grep acen cytoBand.txt | awk ' NR % 2 == 0' | sed 's/chr//' | cut -f 1-2 > ucsc_centromere_mid_hg19.txt
centromere_positions = read.table(CENTROMERE_FILE, header=F, stringsAsFactors=F)
colnames(centromere_positions) = c("chromosome", "position")

merge_complete = function(cndata, levels_remove, allowed_gap, centromere_positions) {
  # Set a flag to denote whether a segment is a composite
  cndata$is_merged = F
  new_cndata = data.frame()
  for (chrom in unique(cndata$chromosome)) {
    print(chrom)
    cndata_chrom = cndata[cndata$chromosome==chrom,]
    
    centromere_boundary = centromere_positions$position[centromere_positions$chromosome==chrom]
    boundary_segment = cndata_chrom$start <= centromere_boundary & cndata_chrom$end >= centromere_boundary
    
    if (sum(boundary_segment)==0 & cndata_chrom$end[nrow(cndata_chrom)] > centromere_boundary) {
      # Take the segment before the centromere boundary
      boundary_segment = which(cndata_chrom$end > centromere_boundary)[1] - 1
    } else if (cndata_chrom$end[nrow(cndata_chrom)] < centromere_boundary) {
      boundary_segment = nrow(cndata_chrom)
    } else {
      boundary_segment = which(boundary_segment)
    }
    
    cndata_first = cndata_chrom[1:boundary_segment,]
    if (boundary_segment < nrow(cndata_chrom)) {
      cndata_second = cndata_chrom[(boundary_segment+1):nrow(cndata_chrom),]
    } else {
      cndata_second = NULL
    }

    res_first = merge_complete_inner(cndata_first, levels_remove, allowed_gap)
    cndata_first_merged = res_first$cndata
    res_second = merge_complete_inner(cndata_second, levels_remove, allowed_gap)
    cndata_second_merged = res_second$cndata

    # compare last-first segment with first-second segment
    merge_boundary_segments = F
    if (!res_first$last_segment_filtered & !res_second$first_segment_filtered) {
      first = paste0(cndata_first_merged$major_cn[nrow(cndata_first_merged)], "_", cndata_first_merged$minor_cn[nrow(cndata_first_merged)])
      second = paste0(cndata_second_merged$major_cn[1], "_", cndata_second_merged$minor_cn[1])
      
      if (first==second) {
        merge_boundary_segments = T
      }
    }
      
    # Concatenate the output for this chromosome
    new_cndata_chrom = NULL
    if (merge_boundary_segments) {
      merged = merge_segments(rbind(cndata_first_merged[nrow(cndata_first_merged),,drop=F], cndata_second_merged[1,,drop=F]), 1, 2)
      
      if (nrow(cndata_first_merged)==1) {
        list_of_dfs = list(merged)
      } else {
        list_of_dfs = list(cndata_first_merged[1:(nrow(cndata_first_merged)-1),,drop=F],
                           merged)
      }
      
      if (nrow(cndata_second_merged) > 1) {
        list_of_dfs[[length(list_of_dfs)+1]] = cndata_second_merged[2:nrow(cndata_second_merged),,drop=F]
      }
      new_cndata_chrom = do.call(rbind, list_of_dfs)
      
      
    } else if (!is.null(cndata_first_merged) & !is.null(cndata_second_merged)) {
      new_cndata_chrom = rbind(cndata_first_merged, cndata_second_merged)
    } else if (!is.null(cndata_first_merged)) {
      new_cndata_chrom = rbind(cndata_first_merged)
    } else if (!is.null(cndata_second_merged)) {
      new_cndata_chrom = rbind(cndata_second_merged)
    } 
    
    # Add to the total segments for this sample
    new_cndata = rbind(new_cndata, new_cndata_chrom)
  }
  return(new_cndata)
}

merge_complete_inner = function(cndata, levels_remove, allowed_gap) {
  cndata_filtered = cndata[cndata$level %in% levels_remove,, drop=F]
  if (is.null(cndata) || nrow(cndata_filtered) == nrow(cndata)) {
	return(list(cndata=NULL, first_segment_filtered=TRUE, last_segment_filtered=TRUE))
  }
  # Keep track of these to provide for a higher-up merging step
  first_segment_filtered = 1 %in% which(cndata$level %in% levels_remove)
  last_segment_filtered = nrow(cndata) %in% which(cndata$level %in% levels_remove)
  cndata = cndata[!cndata$level %in% levels_remove,, drop=F]
  
  merged = T
  while (merged) {
    # If cndata consists of a single segment, then do not proceed
    if (nrow(cndata) <= 1) { break }
    
    merged = F
    new_cndata = data.frame()
    for (i in 1:(nrow(cndata)-1)) {
      ############################################################
      # Check of this and next segment are equal
      ############################################################
      pair_equal = paste0(cndata$major_cn[i], "_", cndata$minor_cn[i]) == paste0(cndata$major_cn[i+1], "_", cndata$minor_cn[i+1])
      
      ############################################################
      # Check for any removed segments
      ############################################################
      i_end = cndata$end[i]
      i_next_start = cndata$start[i+1]
      any_filtered = cndata_filtered$start > i_end & cndata_filtered$end < i_next_start
      if (any(any_filtered)) {
        # Check whether the segment(s) in between have the same major/minor states
        filtered_cn = paste0(cndata_filtered$major_cn[any_filtered], "_", cndata_filtered$minor_cn[any_filtered])
        i_equal = all(filtered_cn==filtered_cn[1]) & paste0(cndata$major_cn[i], "_", cndata$minor_cn[i])==filtered_cn[1]
        i_next_equal = all(filtered_cn==filtered_cn[1]) & paste0(cndata$major_cn[i+1], "_", cndata$minor_cn[i+1])==filtered_cn[1]
      } else {
        filtered_cn = NA
        i_equal = F
        i_next_equal = F
      }
      
      ############################################################
      # Logic to determine whether merging is in order
      ############################################################
      # Check for a large gap between the segments
      large_gap = (i_next_start - i_end) > allowed_gap

      if (!pair_equal) { merge_pair = F }
      
      if (pair_equal & !large_gap & !any(any_filtered)) { merge_pair = T }
            
      if (pair_equal & large_gap & !any(any_filtered)) { merge_pair = F }
      
      if (pair_equal & any(any_filtered) & (!i_equal | !i_next_equal)) { merge_pair = F }
      
      if (pair_equal & any(any_filtered) & (i_equal & i_next_equal)) { merge_pair = T }
      
      ############################################################
      # Perform actions depending on the logic
      ############################################################
      if (merge_pair) {
        new = merge_segments(cndata, i, i+1)
        merged = T
        
        # Append the new segment plus anything that was not assessed still
        new_cndata = rbind(new_cndata, new)
        if (nrow(cndata) > i+1) {
          new_cndata = rbind(new_cndata, cndata[(i+2):nrow(cndata),,drop=F])
        }
        
        # break the for loop at this stage and start over
        cndata = new_cndata
        break
      } else {
        new_cndata = rbind(new_cndata, cndata[i,,drop=F])
        # If this is the final iteration of the for loop and we haven't merged, then add the final i+1 segment as well
        if (i==(nrow(cndata)-1)) {
          new_cndata = rbind(new_cndata, cndata[i+1,,drop=F])
        }
      }
    }
    cndata = new_cndata
  }
  
  #' new_cndata contains all the new segments and should be returned here, not cndata because that is only 
  #' updated if a pair of segments is merged. The final iteration of the for loop therefore appends all 
  #' final segments to new_cndata
  return(list(cndata=cndata, first_segment_filtered=first_segment_filtered, last_segment_filtered=last_segment_filtered))
}

merge_segments = function(cndata, i, j) {
  new = cndata[i,,drop=F]
  new$end = cndata$end[j]
  new$is_merged = T
  return(new)
}

test_merge = function() {
  # set this to 2 for not WGD and 4 to WGD samples
  ploidy = 2
  is_male = TRUE
  # bb_subclones_file = "battenberg/0040b1b6-b07a-4b6e-90ef-133523eaf412_subclones.txt"
  consensus_file = "battenberg/0040b1b6-b07a-4b6e-90ef-133523eaf412.consensus.20170119.somatic.cna.annotated.txt"
  samplename = "0040b1b6-b07a-4b6e-90ef-133523eaf412"
  
  # get centromere mid: grep acen cytoBand.txt | awk ' NR % 2 == 0' | sed 's/chr//' | cut -f 1-2 > ucsc_centromere_mid_hg19.txt
  centromere_positions = read.table("~/repo/icgc_paper_figures/ucsc_centromere_mid_hg19.txt", header=F, stringsAsFactors=F)
  colnames(centromere_positions) = c("chromosome", "position")
  
  levels_remove = c("g", "h")
  allowed_gap = 500000
  
  # Read in the data
  # cndata = read.table(bb_subclones_file, header=T, stringsAsFactors=F)
  cndata_raw = read.table(consensus_file, header=T, stringsAsFactors=F)
  
  #######################################################
  # Test merging of segments and merging of chrom arms
  #######################################################
  cndata = cndata_raw[cndata_raw$chromosome=="1",]
  cndata_temp = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
  if (nrow(cndata_temp) != 2) { print("Test 1 failed") }
  
  #######################################################
  # Test filter with equal CN state segment removed
  #######################################################
  cndata = cndata_raw[cndata_raw$chromosome=="1",]
  cndata$level[7] = "h"
  cndata_temp = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
  if (nrow(cndata_temp) != 2) { print("Test 2 failed") }
  
  #######################################################
  # Test filter with different CN state segment removed
  #######################################################
  cndata = cndata_raw[cndata_raw$chromosome=="1",]
  cndata$level[7] = "h"
  cndata$major_cn[7] = 3
  cndata$total_cn[7] = cndata$major_cn[7] + cndata$minor_cn[7]
  cndata_temp = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
  if (nrow(cndata_temp) != 3) { print("Test 3 failed") }
  
  #######################################################
  # Test large gap
  #######################################################
  cndata = cndata_raw[cndata_raw$chromosome=="1",]
  cndata = cndata[-7,]
  cndata_temp = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
  if (nrow(cndata_temp) != 3) { print("Test 4 failed") }
  
  #######################################################
  # Test filter centromere
  #######################################################
  cndata = cndata_raw[cndata_raw$chromosome=="1",]
  cndata$level[4] = "h"
  cndata_temp = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
  if (nrow(cndata_temp) != 3) { print("Test 4 failed") }
}



merge_simple = function(cndata) {
  new_cndata = data.frame()
  for (chrom in unique(cndata$chromosome)) {
    cndata_chrom = cndata[cndata$chromosome==chrom,]
    inven = paste0(cndata_chrom$major_cn, "_", cndata_chrom$minor_cn)
    a = rle(inven)
    
    start = 0 # start always contains the end position of the previous segment
    for (i in 1:length(a$lengths)) {
      end = start+a$lengths[i]
      new_cndata = rbind(new_cndata, data.frame(chromosome=chrom, start=cndata_chrom$start[start+1], end=cndata_chrom$end[end], major_cn=cndata_chrom$major_cn[end], minor_cn=cndata_chrom$minor_cn[end], frac1_A=1, nMaj2_A=NA, nMin2_A=NA, frac2_A=NA))
      start = end
    }
  }
  return(new_cndata)
}


classify_segments = function(cndata, is_male, ploidy, samplename) {
  
  if (nrow(cndata)==0) {
    return(NULL)
  }
  
  if (is_male) {
    x_exp = ploidy/2; y_exp = ploidy/2; exp_sex_chrom_lvl = ploidy/2;
  } else {
    x_exp = ploidy; y_exp = 0; exp_sex_chrom_lvl = 0;
  }
  
  # If there is no column with a tumour name, add it in temporarily
  if (! "Tumour_Name" %in% colnames(cndata)) {
    cndata = data.frame(Tumour_Name=samplename, cndata)
  }
  
  # allsegs = data.frame(cndata[, c("Tumour_Name", "chromosome", "start", 
  #                                 "end", "major_cn", "minor_cn", 
  #                                 "frac1_A", "nMaj2_A", "nMin2_A", 
  #                                 "frac2_A", "SDfrac_A")], tumour_ploidy=ploidy)
  allsegs = data.frame(cndata[, c("Tumour_Name", "chromosome", "start", 
                                  "end", "major_cn", "minor_cn")], tumour_ploidy=ploidy)
  
  # Now classify all segments into a category
  tot = nrow(allsegs)
  # if you have clonal LOH, subclonal LOH is not counted!!!  Losses not counted
  allsegsa <- NULL
  CNA <- NULL
  for (i in 1:dim(allsegs)[1]) {
	  print(i)
	  save(file="test.RData", i, allsegs)
    select_columns = c(1:7)
    
    # Save some states
    is_hd = allsegs$minor_cn[i] == 0 & allsegs$major_cn[i] == 0  
    is_loh = xor(allsegs$minor_cn[i] == 0, allsegs$major_cn[i] == 0)
    
    
    # These are the same for both types of chromosomes
    # HD
    if (is_hd) {
      allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
      CNA <- c(CNA, "cHD")
    }
    
    # LOH
    if (is_loh) {
      allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
      CNA <- c(CNA, "cLOH")
    }
    
    # Clonal copy number
    if (is_male & allsegs$chromosome[i] %in% c("X", "Y")) {
      
      ####################################################    
      # Sex Chromosomes
      ####################################################
      
      # Gains - male
      # TODO is it strictly required here to check for the other allele to be 0, can this be done upstream?
      if ((allsegs$major_cn[i] > exp_sex_chrom_lvl) | (allsegs$minor_cn[i] > exp_sex_chrom_lvl)) {
        # if ((((allsegs$major_cn[i] > exp_sex_chrom_lvl*3) & (allsegs$minor_cn[i] == 0)) | ((allsegs$major_cn[i] == 0) & (allsegs$minor_cn[i] > exp_sex_chrom_lvl*3)))) {
        #   allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        #   CNA <- c(CNA, "cAmp")
        # }
        # else {
          allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
          CNA <- c(CNA, "cGain")
        # }
      }
      
      
      # Not aberrated - male
      if ( ((allsegs$major_cn[i] == exp_sex_chrom_lvl) & (allsegs$minor_cn[i] == 0)) | ((allsegs$major_cn[i] == 0) & (allsegs$minor_cn[i] == exp_sex_chrom_lvl)) ) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        CNA <- c(CNA, "NoCNA")
      }
      
      
      # Loss - male
      if (((allsegs$major_cn[i] < exp_sex_chrom_lvl) & (allsegs$minor_cn[i] == 0)) | ((allsegs$minor_cn[i] < exp_sex_chrom_lvl) & (allsegs$major_cn[i] == 0))) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        CNA <- c(CNA, "cLoss")
      }
      
    } else {
      ####################################################    
      # Autosomes
      ####################################################
      
      # Gains
      if ((allsegs$major_cn[i] > ploidy/2) | (allsegs$minor_cn[i] > ploidy/2)) {
        # if ((allsegs$major_cn[i] > ploidy*2) | (allsegs$minor_cn[i] > ploidy*2)) {
        #   allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        #   CNA <- c(CNA, "cAmp")
        # }
        # else {
          allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
          CNA <- c(CNA, "cGain")
        # }
      }
      
      
      # Not aberrated
      if ((allsegs$major_cn[i] == ploidy/2) & (allsegs$minor_cn[i] == ploidy/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        CNA <- c(CNA, "NoCNA")
      }
      
      
      # Loss
      if ((allsegs$major_cn[i] < ploidy/2) | (allsegs$minor_cn[i] < ploidy/2)) {
        allsegsa <- rbind(allsegsa, allsegs[i,select_columns])
        CNA <- c(CNA, "cLoss")
      }
    }
  }
    
    
  allsegsa <- cbind(allsegsa, CNA)
  allsegsa$CNA = factor(allsegsa$CNA, levels=c("cGain", "cLoss", "cHD", "NoCNA", "cLOH", "WGD"))
  # allsegsa$frac1_A[allsegsa$CNA == "cGain"|allsegsa$CNA == "cAmp"|allsegsa$CNA == "cLOH"|allsegsa$CNA == "cHD"|allsegsa$CNA == "cLoss"] = 1
  
  # Switch the columns
  allsegsa = data.frame(allsegsa[,-1], Tumour_Name=samplename)
  
  return(allsegsa)
}

# Read in the data
# cndata = read.table(bb_subclones_file, header=T, stringsAsFactors=F)
cndata = read.table(consensus_file, header=T, stringsAsFactors=F)

# Remove all segments witout a call
if (any(is.na(cndata$minor_cn) | is.na(cndata$major_cn))) {
  cndata = cndata[!(is.na(cndata$minor_cn) | is.na(cndata$major_cn)),]
}

cndata_simple = merge_simple(cndata)
cndata_simple = classify_segments(cndata_simple, is_male, ploidy, samplename)

levels_remove = c("g", "h")
allowed_gap = 500000
cndata_complete = merge_complete(cndata, levels_remove, allowed_gap, centromere_positions)
cndata_complete = classify_segments(cndata_complete, is_male, ploidy, samplename)

if (!is.null(cndata_simple)) {
  summ_simple = table(cndata_simple$CNA)
  a = data.frame(t(matrix(as.numeric(summ_simple))))
  colnames(a) = paste0("simple_", names(summ_simple))
} else {
  a = data.frame(cGain=0, cLoss=0, cHD=0, NoCNA=0, cLOH=0, WGD=0)
  colnames(a) = paste0("simple_", colnames(a))
}

if (!is.null(cndata_complete)) {
  summ_complete = table(cndata_complete$CNA)
  b = data.frame(t(matrix(as.numeric(summ_complete))))
  colnames(b) = paste0("complete_", names(summ_complete))
} else {
  b = data.frame(cGain=0, cLoss=0, cHD=0, NoCNA=0, cLOH=0, WGD=0)
  colnames(b) = paste0("simple_", colnames(b))
}

output = data.frame(samplename=samplename, a, simple_wgd=ifelse(is_wgd, 1, 0), b, complete_wgd=ifelse(is_wgd, 1, 0))

write.table(output, file=file.path(outdir, paste0(samplename, "_cna_eventcount.txt")), sep="\t", quote=F, row.names=F)
save.image(file=file.path(outdir, paste0(samplename, "_cna_eventcount.RData")))

