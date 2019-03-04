
# load in CLIPdb for human and mouse
# what other RBPs bind the FUS retained introns other than FUS itself?
library(dplyr)
library(readr)
library(ggplot2)

## ---------- mouse ------------------------------------------------------------

mouse_clipdb <- 
  read_tsv("~/Downloads/mouse_RBP_binding_sites.txt", 
           col_names = c("chr", "start", "end", "info", "V6", 
                         "strand", "RBP", "protocol",
                         "tissue", "accession", "score"
           )
  )

mouse_clipdb$length <- mouse_clipdb$end - mouse_clipdb$start

fus_coords_mouse <- list(chr = "chr7", start = 127965479, end = 127984031)
fus_coords_mouse_intron <- list(chr = "chr7", start = 127972766, end = 127975906)

fus_mouse <- 
  filter(mouse_clipdb, 
         chr == fus_coords_mouse_intron$chr, 
         start > fus_coords_mouse_intron$start, 
         end < fus_coords_mouse_intron$end
         )

# compute per-base coverage for each RBP

mouse_fus_coverage <- 
  group_by(fus_mouse, RBP, protocol, tissue, accession) %>%
  summarise( fus_n = n(), fus_coverage = sum(length) ) 

# background: coverage over the entire mouse genome
mouse_background_coverage <-
  group_by(mouse_clipdb, RBP, protocol, tissue, accession) %>%
  summarise( background_n = n(), background_coverage = sum(length) ) 

mouse_combined <- left_join(mouse_fus_coverage, mouse_background_coverage)

mouse_combined$cov_enrichment <- mouse_combined$fus_coverage / mouse_combined$background_coverage

ggplot(mouse_combined, aes(x = RBP, y = cov_enrichment) ) + geom_col( aes(fill = protocol)) +
  coord_flip() +
  facet_wrap(~protocol, scales = "free")

#cor.test(combined$fus_coverage, combined$background_coverage)

ggplot(mouse_combined, aes(x = log(background_coverage), y = log(fus_coverage), label = RBP)) +
  geom_text(aes(colour = protocol))

ggplot(mouse_combined, aes(x = log(background_n), y = log(fus_n), label = RBP)) +
  geom_text(aes(colour = protocol))


filter(mouse_combined, RBP %in% c("TARDBP", "FUS" ,"TAF15", "EWSR1", "APC", "U2AF2", "SRSF1")) %>%
  split( 1:nrow(.) )  %>%
  purrr::map( ~{
    
    outfile <-  paste0("data/iCLIP/usual_suspects/mouse/", .x$RBP, "_", gsub(",", "_",.x$protocol ), "_",.x$tissue, ".bedGraph" )
    print(outfile)
    peaks <- filter(fus_mouse, RBP == .x$RBP, protocol == .x$protocol, tissue == .x$tissue) %>%
      select( chr, start, end, score)
    write_tsv(peaks, path = outfile, col_names = FALSE)
  })

# APC replicates

filter(mouse_combined, RBP == c("APC")) %>%
  split( 1:nrow(.) )  %>%
  purrr::map( ~{
    
    outfile <-  paste0("data/iCLIP/usual_suspects/mouse/",
                       .x$RBP, "_", 
                       gsub(",", "_",.x$protocol ), "_",
                       .x$tissue, "_",
                       .x$accession, ".bedGraph" )
    print(outfile)
    peaks <- filter(fus_mouse, RBP == .x$RBP, protocol == .x$protocol, tissue == .x$tissue, accession == .x$accession) %>%
      select( chr, start, end, score)
    write_tsv(peaks, path = outfile, col_names = FALSE)
  })

# now check coverage over set of retained introns and compare to a null.
load("data/SGSeq/joint_analysis/novel_events/processed_all_joint.Rdata")
introns <- filter(res.clean, status_relaxed == "Overlapping", type == "RI")
nulls <- filter(res.clean, status_relaxed == "None", type == "RI")

# make into GRange
mouse_clipdb_grange <- GenomicRanges::GRanges(
  seqnames = mouse_clipdb$chr,
  ranges = IRanges(start = mouse_clipdb$start,
                   end = mouse_clipdb$end),
  RBP = mouse_clipdb$RBP,
  protocol = mouse_clipdb$protocol,
  tissue = mouse_clipdb$tissue,
  accession = mouse_clipdb$accession,
  length = mouse_clipdb$length
)

# find coverage over set of ranges
clip_coverage <- function(coord){
  # split coords
  coords <- str_split_fixed(coord, pattern = ":|-", n = 3)[1,]

  
  # get clipdb entries within that range
  range <-
   filter(mouse_clipdb,
          chr == coords[1],
          start > coords[2],
          end < coords[3]
   )
    # compute per-base coverage for each RBP
  
  coverage <- 
    group_by(range, RBP, protocol, tissue, accession) %>%
    summarise( n = n(), coverage = sum(length) ) 
  
  return(coverage)
}

coord <- introns$coords[1:10]

clip_coverage_grange <- function(input_coords){
  # split coords
  coords <- str_split_fixed(input_coords, pattern = ":|-", n = 3)
  
  coords <- list(chr = coords[,1],
                 start = as.numeric(coords[,2]),
                 end = as.numeric(coords[,3]))
  
  coords_range <- GRanges(seqnames = coords$chr, 
                          ranges = IRanges(as.numeric(coords$start),
                                           as.numeric(coords$end) )
  )
  
  overlap <- findOverlaps(mouse_clipdb_grange, coords_range)
  range <- mouse_clipdb_grange[ queryHits(overlap)] %>% as.data.frame()
  
  
  # split by input intron
  range_split <- split(range, subjectHits(overlap))
  
  # nest dataframe
  coverage_list <- 
    range %>%
    mutate( intron = subjectHits(overlap)) %>%
    group_by(intron) %>%
    nest() %>%
    mutate( summary = map(
      data, ~{
        .x %>%
          group_by(RBP, protocol, tissue, accession) %>%
          summarise( n = n(), coverage = sum(length) ) 
      }
    ))
  
  # match back original coords
  intron_coords <- tibble(
    intron = 1:length(input_coords),
    coords = input_coords
  )
  
  coverage_list <- coverage_list %>%
    full_join(intron_coords, by = "intron")
  
  # compute per-base coverage for each RBP for each intron


  return(coverage_list)
#  return(range)  
}


introns_coverage <- clip_coverage_grange(introns$coords)

null_coverage <- clip_coverage_grange(nulls$coords[1:1000])

# for each intron, how many have FUS binding?
map_lgl(introns_coverage$summary, ~{
  if( is.null(.x) ){return(FALSE)}
  "FUS" %in% .x$RBP
}) %>%
  table()

# do for all RBPs in data
all_rbps <- mouse_clipdb %>%
  select(RBP, tissue, protocol, accession) %>%
  distinct() %>%
  split(1:nrow(.))

# for each RBP go through each set
results <- 
  all_rbps %>%
  map_df(~{
    df <- .x
    # iterate through 
    x_introns <- 
      map_lgl(introns_coverage$summary, ~{
        if( is.null(.x) ){return(FALSE)}
        # find whether each RBP is  present in any intron
        if_true <- filter(.x, RBP == df$RBP, tissue == df$tissue, protocol == df$protocol, accession == df$accession)
        return( nrow(if_true) > 0 )
    }) %>% sum()
    
    x_nulls <- 
      map_lgl(null_coverage$summary, ~{
        if( is.null(.x) ){return(FALSE)}
        # find whether each RBP is  present in any intron
        if_true <- filter(.x, RBP == df$RBP, tissue == df$tissue, protocol == df$protocol, accession == df$accession)
        return( nrow(if_true) > 0 )
      }) %>% sum()
    
    df %>%
      mutate( x_introns, x_nulls)
  })



## ---------- human


# human CLIPdb
human_clipdb <- read_tsv("~/Downloads/human_RBP_binding_sites.txt",
                         col_names = c("chr", "start", "end", "info", "V6", 
                                       "strand", "RBP", "protocol",
                                       "tissue", "accession", "score"
                         )         
)

human_clipdb$length <- human_clipdb$end - human_clipdb$start


fus_coords_human_intron <- list(chr = "chr16", 
                                start = 31185155, 
                                end = 31188343)


fus_human <- 
  filter(human_clipdb, 
         chr == fus_coords_human_intron$chr, 
         start > fus_coords_human_intron$start, 
         end < fus_coords_human_intron$end
  )

fus_coverage <- 
  group_by(fus_human, RBP, protocol, tissue, accession) %>%
  summarise( fus_n = n(), fus_coverage = sum(length) ) 

background_coverage <-
  group_by(human_clipdb, RBP, protocol, tissue, accession) %>%
  summarise( background_n = n(), background_coverage = sum(length) ) 

combined <- left_join(fus_coverage, background_coverage)

combined$cov_enrichment <- combined$fus_coverage / combined$background_coverage

# plot coverage of an RBP over the FUS intron
fus_intron_length <- fus_coords_human_intron$end - fus_coords_human_intron$start

combined$fus_coverage_prop <- combined$fus_coverage / fus_intron_length

ggplot(combined, aes( x = log(background_n), y = fus_n)) +
  geom_text( aes(label = RBP, colour = protocol))

# cleavage factors are present in the human FUS intron
# where do they bind? Export out peaks to bed file
cleavage_peaks <- filter(fus_human, RBP == "CPSF6" | RBP == "CSTF2")
# get just iclip

cleavage_peaks <- filter(cleavage_peaks, grepl("iCLIP,CIMS", protocol)) %>%
  select(chr, start, end, score)
write_tsv(cleavage_peaks, path = "data/iCLIP/cleavage_factors_human_clip.bedGraph",col_names=FALSE)

# for the top 20 factors - write out bedGraph files
arrange(combined, desc(cov_enrichment) ) %>% head(50) %>%
  split( 1:nrow(.) )  %>%
  purrr::map( ~{

    outfile <-  paste0("data/iCLIP/", .x$RBP, "_", gsub(",", "_",.x$protocol ), "_",.x$tissue, ".bedGraph" )
    print(outfile)
    peaks <- filter(fus_human, RBP == .x$RBP, protocol == .x$protocol, tissue == .x$tissue) %>%
      select( chr, start, end, score)
    write_tsv(peaks, path = outfile, col_names = FALSE)
  })

filter(combined, RBP %in% c("TARDBP", "FUS" ,"TAF15", "EWSR1")) %>%
  split( 1:nrow(.) )  %>%
  purrr::map( ~{
    
    outfile <-  paste0("data/iCLIP/usual_suspects/", .x$RBP, "_", gsub(",", "_",.x$protocol ), "_",.x$tissue, ".bedGraph" )
    print(outfile)
    peaks <- filter(fus_human, RBP == .x$RBP, protocol == .x$protocol, tissue == .x$tissue) %>%
      select( chr, start, end, score)
    write_tsv(peaks, path = outfile, col_names = FALSE)
  })

