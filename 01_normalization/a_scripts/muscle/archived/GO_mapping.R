config <- list(
  # Data files
  go_map_muscle = "00_input/muscle/CRm_raw_map.csv",
  # Directories 
  file_dir = "00_input/muscle"
)
# Data restructuring and mapping

raw_map <- read.csv(file.path(config$go_map_muscle))
head(raw_map)

raw_map <- raw_map %>%
  # Split by ";" into a list column
  mutate(comp_list = str_split(Gene.Ontology..cellular.component., ";")) %>%
  # Remove GO IDs and trim whitespace
  mutate(comp_list = lapply(comp_list, function(x) str_trim(str_remove(x, "\\s*\\[GO:\\d+\\]"))))

raw_map <- raw_map %>%
  mutate(
    comp1 = sapply(comp_list, function(x) if(length(x) >= 1) x[1] else NA),
    comp2 = sapply(comp_list, function(x) if(length(x) >= 2) x[2] else NA),
    comp3 = sapply(comp_list, function(x) if(length(x) >= 3) x[3] else NA)
  )

comp_map <- tibble(
  pattern = c("acrosomal", 
              "actin", 
              "axon", 
              "autolysosome",
              "autophagosome", 
              "azurophil",
              "blood microparticle",
              "brush border",
              "calcium",
              "cardiac",
              "chromosome",
              "Golgi",
              "clathrin",
              "collagen",
              "COP",
              "Cul",
              "cytoplasm",
              "cytosol",
              "endosome",
              "endoplasmic reticulum",
              "extracellular",
              "lysosomal",
              "mitochondrial",
              "nuclear",
              "sarcoplasmic reticulum"),
  category = c("acrosomal", 
               "actin", 
               "axon", 
               "autolysosome",
               "autophagosome", 
               "azurophil",
               "blood microparticle",
               "brush border",
               "calcium",
               "cardiac",
               "chromosome",
               "Golgi",
               "clathrin",
               "collagen",
               "COP",
               "Cul",
               "cytoplasm",
               "cytosol",
               "endosome",
               "endoplasmic reticulum",
               "extracellular",
               "lysosomal",
               "mitochondrial",
               "nuclear",
               "sarcoplasmic reticulum")
)

condense_compartment <- function(x) {
  dplyr::case_when(
    # Cytoplasm / Cytosol
    grepl("cytoplasm|cytosol", x, ignore.case = TRUE) ~ "Cytoplasm/Cytosol",
    
    # Mitochondria
    grepl("mitochondr|MICOS complex|m-AAA|endopeptidase Clp", x, ignore.case = TRUE) ~ "Mitochondria",
    
    # Skeletal muscle-specific
    grepl("sarcomere|sarcoplasmic reticulum|contractile muscle fiber|cardiac myofibril|M band|A band|costamere|longitudinal sarcoplasmic|cardiac Troponin|ankyrin-1|dystrophin", x, ignore.case = TRUE) ~ "Sarcoplasmic/Myofibrillar",  # + dystrophin
    
    # Actin / Cytoskeleton
    grepl("actin|cytoskeleton|actomyosin|intermediate filament|microtubule", x, ignore.case = TRUE) ~ "Cytoskeleton/Contractile",  # + intermediate filament, microtubule
    
    # Blood microparticle
    grepl("blood microparticle", x, ignore.case = TRUE) ~ "Blood Microparticle",
    
    # Extracellular / Secreted
    grepl("extracellular|basement membrane|collagen|elastic fiber|discoidal|extracellular region|extracellular space", x, ignore.case = TRUE) ~ "Extracellular/Secreted",  # + elastic fiber
    
    # Membrane / Cell surface
    grepl("membrane|cell surface|caveola|plasma membrane|cell cortex|cell periphery|cell leading edge", x, ignore.case = TRUE) ~ "Membrane/Cell Surface",
    
    # Nucleus / Chromatin
    grepl("nucleus|chromatin|chromosome|nucleoplasm|Cajal|fibrillar|Barr body|CENP|nuclear speck|cleavage body|chromatoid body", x, ignore.case = TRUE) ~ "Nucleus/Chromatin",  # + nuclear speck, cleavage body, chromatoid body
    
    # ER / Golgi / Vesicle
    grepl("endoplasmic reticulum|COPII|COPI|Golgi|endomembrane|annulate|endosome|clathrin|amphisome|autophagosome|autolysosome|lysosom|endocytic|early endosome", x, ignore.case = TRUE) ~ "ER/Golgi/Vesicle",  # + COPI, Golgi apparatus, endosome, clathrin, amphisome, autophagosome variants
    
    # Junctions / Adhesion
    grepl("junction|adherens|anchoring|focal adhesion|tight junction|cell-cell|fascia adherens|contact zone", x, ignore.case = TRUE) ~ "Cell Junction/Adhesion",  # + fascia adherens, contact zone
    
    # Granules / Lysosomes
    grepl("granule|aggresome|vesicle|lamellar body|stress granule", x, ignore.case = TRUE) ~ "Granule/Lysosome/Vesicle",  # + stress granule; lysosome/autophagosome moved to ER/Golgi/Vesicle above
    
    # Cilia / Centrosome / Neuron
    grepl("cili|centrosome|centriole|axoneme|basal body|axolemma|axon|dendrite|synapse|GABA|cerebellar|postsynaptic|cleavage furrow", x, ignore.case = TRUE) ~ "Cilia/Centrosome/Neuron",  # + postsynaptic, cleavage furrow
    
    # Named complexes
    grepl("complex|spliceosome|signalosome|ubiquitin|ligase|chaperonin|ATPase|synthase|synthetase|isocitrate|flavoprotein|glycine cleavage|calcineurin|calpain|CD40|angiogenin|AnxA2|MIB|CRD|iron-sulfur|aminoacyl|calmodulin|kinase|Arp2|beta-catenin|Cul[0-9]|aryl hydrocarbon|AP-[0-9]|BRISC|BRCA|branched-chain|acetolactate|chloride channel|commitment|core-binding|Cdc73|DNA replication|snoRNP|bub[0-9]|CAAX|B-WICH|Atg|adaptor", x, ignore.case = TRUE) ~ "Protein Complex",  # + AP-1/2, BRISC, BRCA1-A, branched-chain dehydrogenase, chloride channel, and other new named complexes
    
    # Acrosome
    grepl("acrosom", x, ignore.case = TRUE) ~ "Acrosome",
    
    # Empty
    x == "" ~ "Unknown",
    
    TRUE ~ "Other"
  )
}

# Build GO mapping (same logic you already used)
go_map <- raw_map %>%
  dplyr::select(From, comp1, comp2, comp3) %>%
  mutate(
    comp1_grouped = condense_compartment(comp1)
  )

write.csv(go_map, file.path(config$file_dir, "GO_map.csv"))
