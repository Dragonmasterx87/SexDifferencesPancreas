# NAalysis of Ruth's data
# Hpap dataset
system.time({
  #user    system elapsed 
  #9294.29 614.52 9908.50 (~2.75hrs)
hpap <- readRDS(file=r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\Raw_hpap_file\hpap_islet_scRNAseq.rds)")
table(hpap@meta.data[["Library"]])

Idents(hpap) <- "Library"
hpap$ancestry <- plyr::mapvalues(
  x= hpap$Library,
  from = c("HPAP-019",
           "HPAP-020",
           "HPAP-021",
           "HPAP-022",
           "HPAP-023",
           "HPAP-024",
           "HPAP-026",
           "HPAP-028",
           "HPAP-029",
           "HPAP-032",
           "HPAP-034",
           "HPAP-035",
           "HPAP-036",
           "HPAP-037",
           "HPAP-038",
           "HPAP-039",
           "HPAP-040",
           "HPAP-042",
           "HPAP-043",
           "HPAP-044",
           "HPAP-045",
           "HPAP-047",
           "HPAP-049",
           "HPAP-050",
           "HPAP-051",
           "HPAP-052",
           "HPAP-053",
           "HPAP-054",
           "HPAP-055",
           "HPAP-056",
           "HPAP-057",
           "HPAP-058",
           "HPAP-059",
           "HPAP-061",
           "HPAP-063",
           "HPAP-064",
           "HPAP-065",
           "HPAP-070",
           "HPAP-071",
           "HPAP-072",
           "HPAP-074",
           "HPAP-075",
           "HPAP-077",
           "HPAP-079",
           "HPAP-080",
           "HPAP-081",
           "HPAP-082",
           "HPAP-083",
           "HPAP-084",
           "HPAP-085",
           "HPAP-087",
           "HPAP-088",
           "HPAP-090",
           "HPAP-091",
           "HPAP-092",
           "HPAP-099",
           "HPAP-100",
           "HPAP-101",
           "HPAP-103",
           "HPAP-104",
           "HPAP-105",
           "HPAP-106",
           "HPAP-107",
           "HPAP-108",
           "HPAP-109"),
  to = c("white", #"HPAP-019",
         "white", #"HPAP-020",
         "white", #"HPAP-021",
         "white", #"HPAP-022",
         "white", #"HPAP-023",
         "white", #"HPAP-024",
         "white", #"HPAP-026",
         "white", #"HPAP-028",
         "white", #"HPAP-029",
         "white", #"HPAP-032",
         "white", #"HPAP-034",
         "white", #"HPAP-035",
         "white", #"HPAP-036",
         "white", #"HPAP-037",
         "white", #"HPAP-038",
         "white", #"HPAP-039",
         "white", #"HPAP-040",
         "white", #"HPAP-042",
         "hispanic", #"HPAP-043",
         "white", #"HPAP-044",
         "white", #"HPAP-045",
         "white", #"HPAP-047",
         "white",  #"HPAP-049",
         "hispanic", #"HPAP-050",
         "black", #"HPAP-051",
         "black", #"HPAP-052",
         "white", #"HPAP-053",
         "white", #"HPAP-054",
         "hispanic", #"HPAP-055",
         "white", #"HPAP-056",
         "white", #"HPAP-057",
         "black", #"HPAP-058",
         "white", #"HPAP-059",
         "black", #"HPAP-061",
         "white", #"HPAP-063",
         "black", #"HPAP-064",
         "white", #"HPAP-065",
         "black", #"HPAP-070",
         "white", #"HPAP-071",
         "hispanic", #"HPAP-072",
         "white", #"HPAP-074",
         "white", #"HPAP-075",
         "white", #"HPAP-077",
         "hispanic", #"HPAP-079",
         "hispanic", #"HPAP-080",
         "white", #"HPAP-081",
         "white", #"HPAP-082",
         "black", #"HPAP-083",
         "white", #"HPAP-084",
         "white", #"HPAP-085",
         "white", #"HPAP-087",
         "white", #"HPAP-088",
         "white", #"HPAP-090",
         "hispanic", #"HPAP-091",
         "hispanic", #"HPAP-092",
         "hispanic",  #"HPAP-099",
         "white", #"HPAP-100",
         "hispanic", #"HPAP-101",
         "white", #"HPAP-103",
         "hispanic", #"HPAP-104",
         "hispanic", #"HPAP-105",
         "white", #"HPAP-106",
         "white", #"HPAP-107",
         "black", #"HPAP-108",
         "hispanic" #"HPAP-109"
  )
)

table(hpap@meta.data[["Library"]])
table(hpap@meta.data[["Sex"]])
table(hpap@meta.data[["ancestry"]])

# Subsetting data
Idents(hpap) <- "Library"
{
  HPAP_019 <- subset(hpap, idents = c("HPAP-019"))
  HPAP_020 <- subset(hpap, idents = c("HPAP-020"))
  HPAP_021 <- subset(hpap, idents = c("HPAP-021"))
  HPAP_022 <- subset(hpap, idents = c("HPAP-022"))
  HPAP_023 <- subset(hpap, idents = c("HPAP-023"))
  HPAP_024 <- subset(hpap, idents = c("HPAP-024"))
  HPAP_026 <- subset(hpap, idents = c("HPAP-026"))
  HPAP_028 <- subset(hpap, idents = c("HPAP-028"))
  HPAP_029 <- subset(hpap, idents = c("HPAP-029"))
  HPAP_032 <- subset(hpap, idents = c("HPAP-032"))
  HPAP_034 <- subset(hpap, idents = c("HPAP-034"))
  HPAP_035 <- subset(hpap, idents = c("HPAP-035"))
  HPAP_036 <- subset(hpap, idents = c("HPAP-036"))
  HPAP_037 <- subset(hpap, idents = c("HPAP-037"))
  HPAP_038 <- subset(hpap, idents = c("HPAP-038"))
  HPAP_039 <- subset(hpap, idents = c("HPAP-039"))
  HPAP_040 <- subset(hpap, idents = c("HPAP-040"))
  HPAP_042 <- subset(hpap, idents = c("HPAP-042"))
  HPAP_043 <- subset(hpap, idents = c("HPAP-043"))
  HPAP_044 <- subset(hpap, idents = c("HPAP-044"))
  HPAP_045 <- subset(hpap, idents = c("HPAP-045"))
  HPAP_047 <- subset(hpap, idents = c("HPAP-047"))
  HPAP_049 <- subset(hpap, idents = c("HPAP-049"))
  HPAP_050 <- subset(hpap, idents = c("HPAP-050"))
  HPAP_051 <- subset(hpap, idents = c("HPAP-051"))
  HPAP_052 <- subset(hpap, idents = c("HPAP-052"))
  HPAP_053 <- subset(hpap, idents = c("HPAP-053"))
  HPAP_054 <- subset(hpap, idents = c("HPAP-054"))
  HPAP_055 <- subset(hpap, idents = c("HPAP-055"))
  HPAP_056 <- subset(hpap, idents = c("HPAP-056"))
  HPAP_057 <- subset(hpap, idents = c("HPAP-057"))
  HPAP_058 <- subset(hpap, idents = c("HPAP-058"))
  HPAP_059 <- subset(hpap, idents = c("HPAP-059"))
  HPAP_061 <- subset(hpap, idents = c("HPAP-061"))
  HPAP_063 <- subset(hpap, idents = c("HPAP-063"))
  HPAP_064 <- subset(hpap, idents = c("HPAP-064"))
  HPAP_065 <- subset(hpap, idents = c("HPAP-065"))
  HPAP_070 <- subset(hpap, idents = c("HPAP-070"))
  HPAP_071 <- subset(hpap, idents = c("HPAP-071"))
  HPAP_072 <- subset(hpap, idents = c("HPAP-072"))
  HPAP_074 <- subset(hpap, idents = c("HPAP-074"))
  HPAP_075 <- subset(hpap, idents = c("HPAP-075"))
  HPAP_077 <- subset(hpap, idents = c("HPAP-077"))
  HPAP_079 <- subset(hpap, idents = c("HPAP-079"))
  HPAP_080 <- subset(hpap, idents = c("HPAP-080"))
  HPAP_081 <- subset(hpap, idents = c("HPAP-081"))
  HPAP_082 <- subset(hpap, idents = c("HPAP-082"))
  HPAP_083 <- subset(hpap, idents = c("HPAP-083"))
  HPAP_084 <- subset(hpap, idents = c("HPAP-084"))
  HPAP_085 <- subset(hpap, idents = c("HPAP-085"))
  HPAP_087 <- subset(hpap, idents = c("HPAP-087"))
  HPAP_088 <- subset(hpap, idents = c("HPAP-088"))
  HPAP_090 <- subset(hpap, idents = c("HPAP-090"))
  HPAP_091 <- subset(hpap, idents = c("HPAP-091"))
  HPAP_092 <- subset(hpap, idents = c("HPAP-092"))
  HPAP_099 <- subset(hpap, idents = c("HPAP-099"))
  HPAP_100 <- subset(hpap, idents = c("HPAP-100"))
  HPAP_101 <- subset(hpap, idents = c("HPAP-101"))
  HPAP_103 <- subset(hpap, idents = c("HPAP-103"))
  HPAP_104 <- subset(hpap, idents = c("HPAP-104"))
  HPAP_105 <- subset(hpap, idents = c("HPAP-105"))
  HPAP_106 <- subset(hpap, idents = c("HPAP-106"))
  HPAP_107 <- subset(hpap, idents = c("HPAP-107"))
  HPAP_108 <- subset(hpap, idents = c("HPAP-108"))
  HPAP_109 <- subset(hpap, idents = c("HPAP-109"))
}

pancreas.list <- list("HPAP_019" = HPAP_019, "HPAP_020" = HPAP_020, "HPAP_021" = HPAP_021, "HPAP-022" = HPAP_022, "HPAP_023" = HPAP_023,
                      "HPAP_024" = HPAP_024, "HPAP_026" = HPAP_026, "HPAP_028" = HPAP_028, "HPAP_029" = HPAP_029, "HPAP_032" = HPAP_032,
                      "HPAP_034" = HPAP_034, "HPAP_035" = HPAP_035, "HPAP_036" = HPAP_036, "HPAP_037" = HPAP_037, "HPAP_038" = HPAP_038,
                      "HPAP_039" = HPAP_039, "HPAP_040" = HPAP_040, "HPAP_042" = HPAP_042, "HPAP_043" = HPAP_043, "HPAP_044" = HPAP_044,
                      "HPAP_045" = HPAP_045, "HPAP_047" = HPAP_047, "HPAP_049" = HPAP_049, "HPAP_050" = HPAP_050, "HPAP_051" = HPAP_051,
                      "HPAP_052" = HPAP_052, "HPAP_053" = HPAP_053, "HPAP_054" = HPAP_054, "HPAP_055" = HPAP_055, "HPAP_056" = HPAP_056,
                      "HPAP_057" = HPAP_057, "HPAP_058" = HPAP_058, "HPAP_059" = HPAP_059, "HPAP_061" = HPAP_061, "HPAP_063" = HPAP_063,
                      "HPAP_064" = HPAP_064, "HPAP_065" = HPAP_065, "HPAP_070" = HPAP_070, "HPAP_071" = HPAP_071, "HPAP_072" = HPAP_072,
                      "HPAP_074" = HPAP_074, "HPAP_075" = HPAP_075, "HPAP_077" = HPAP_077, "HPAP_079" = HPAP_079, "HPAP_080" = HPAP_080,
                      "HPAP_081" = HPAP_081, "HPAP_082" = HPAP_082, "HPAP_083" = HPAP_083, "HPAP_084" = HPAP_084, "HPAP_085" = HPAP_085,
                      "HPAP_087" = HPAP_087, "HPAP_088" = HPAP_088, "HPAP_090" = HPAP_090, "HPAP_091" = HPAP_091, "HPAP_092" = HPAP_092,
                      "HPAP_099" = HPAP_099, "HPAP_100" = HPAP_100, "HPAP_101" = HPAP_101, "HPAP_103" = HPAP_103, "HPAP_104" = HPAP_104,
                      "HPAP_105" = HPAP_105, "HPAP_106" = HPAP_106, "HPAP_107" = HPAP_107, "HPAP_108" = HPAP_108, "HPAP_109" = HPAP_109
)

# Pulling data for analysis
pancreas_subset <- pancreas.list[c("HPAP-022", "HPAP_026", "HPAP_035", "HPAP_036", "HPAP_037", "HPAP_040", "HPAP_044", "HPAP_051", 
                                     "HPAP_052", "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", 
                                     "HPAP_063", "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", 
                                     "HPAP_081", "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099",
                                     "HPAP_100", "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109")
]

pancreas_subset <- lapply(X = pancreas_subset, 
                            FUN = function(x){
                              x[['percent.mt']] <- PercentageFeatureSet(x, 
                                                                        pattern = '^MT-')
                              x <- SCTransform(x,
                                               vars.to.regress = "percent.mt", 
                                               verbose = TRUE, 
                                               return.only.var.genes = FALSE, 
                                               vst.flavor = "v2")
                            })


#Save point
#saveRDS(pancreas_subset, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_subset.rds)")

{
  HP2022801 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2022801.rds)")
  SAMN15877725 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\SAMN15877725.rds)")
  HP2024001 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2024001.rds)")
  HP2031401 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2031401.rds)")
  HP2105501 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2105501.rds)")
  HP2106201 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2106201.rds)")
  HP2107001 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2107001.rds)")
  HP2107901 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2107901.rds)")
  HP2108601 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2108601.rds)")
  HP2108901 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2108901.rds)")
  HP2110001 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2110001.rds)")
  HP2121601 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2121601.rds)")
  HP2123201 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2123201.rds)")
  HP2132801 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2132801.rds)")
  HP2202101 <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\CoreArch_CurrentR\HP files after doubletfinder\HP2202101.rds)")
}


# Adjust and add metadata
# Sex specific Metadata addition
{
  HP2022801$Sex <- "M"
  SAMN15877725$Sex <- "M"
  HP2024001$Sex <- "F"
  HP2031401$Sex <- "M"
  HP2105501$Sex <- "F"
  HP2106201$Sex <- "F"
  HP2107001$Sex <- "M"
  HP2107901$Sex <- "M"
  HP2108601$Sex <- "F"
  HP2108901$Sex <- "F"
  HP2110001$Sex <- "M"
  HP2121601$Sex <- "F"
  HP2123201$Sex <- "M"
  HP2132801$Sex <- "F"
  HP2202101$Sex <- "F"
  
  
  HP2022801@meta.data[["Tissue Source"]] <- "Tulane"
  SAMN15877725@meta.data[["Tissue Source"]] <- "Tulane"
  HP2024001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2031401@meta.data[["Tissue Source"]] <- "Tulane"
  HP2105501@meta.data[["Tissue Source"]] <- "Tulane"
  HP2106201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2107901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2108901@meta.data[["Tissue Source"]] <- "Tulane"
  HP2110001@meta.data[["Tissue Source"]] <- "Tulane"
  HP2121601@meta.data[["Tissue Source"]] <- "Tulane"
  HP2123201@meta.data[["Tissue Source"]] <- "Tulane"
  HP2132801@meta.data[["Tissue Source"]] <- "Tulane"
  HP2202101@meta.data[["Tissue Source"]] <- "Tulane"
  
  
  HP2022801@meta.data[["Library"]] <- deparse(substitute(HP2022801))
  SAMN15877725@meta.data[["Library"]] <- deparse(substitute(SAMN15877725))
  HP2024001@meta.data[["Library"]] <- deparse(substitute(HP2024001))
  HP2031401@meta.data[["Library"]] <- deparse(substitute(HP2031401))
  HP2105501@meta.data[["Library"]] <- deparse(substitute(HP2105501))
  HP2106201@meta.data[["Library"]] <- deparse(substitute(HP2106201))
  HP2107001@meta.data[["Library"]] <- deparse(substitute(HP2107001))
  HP2107901@meta.data[["Library"]] <- deparse(substitute(HP2107901))
  HP2108601@meta.data[["Library"]] <- deparse(substitute(HP2108601))
  HP2108901@meta.data[["Library"]] <- deparse(substitute(HP2108901))
  HP2110001@meta.data[["Library"]] <- deparse(substitute(HP2110001))
  HP2121601@meta.data[["Library"]] <- deparse(substitute(HP2121601))
  HP2123201@meta.data[["Library"]] <- deparse(substitute(HP2123201))
  HP2132801@meta.data[["Library"]] <- deparse(substitute(HP2132801))
  HP2202101@meta.data[["Library"]] <- deparse(substitute(HP2202101))
  
  
  HP2022801@meta.data[["sample"]] <- NULL
  SAMN15877725@meta.data[["sample"]] <- NULL
  HP2024001@meta.data[["sample"]] <- NULL
  HP2031401@meta.data[["sample"]] <- NULL
  HP2105501@meta.data[["sample"]] <- NULL
  HP2106201@meta.data[["sample"]] <- NULL
  HP2107001@meta.data[["sample"]] <- NULL
  HP2107901@meta.data[["sample"]] <- NULL
  HP2108601@meta.data[["sample"]] <- NULL
  HP2108901@meta.data[["sample"]] <- NULL
  HP2110001@meta.data[["sample"]] <- NULL
  HP2121601@meta.data[["sample"]] <- NULL
  HP2123201@meta.data[["sample"]] <- NULL
  HP2132801@meta.data[["sample"]] <- NULL
  HP2202101@meta.data[["sample"]] <- NULL
  
  
  HP2022801@meta.data[["Chemistry"]] <- "10Xv3"
  SAMN15877725@meta.data[["Chemistry"]] <- "10Xv3"
  HP2024001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2031401@meta.data[["Chemistry"]] <- "10Xv3"
  HP2105501@meta.data[["Chemistry"]] <- "10Xv3"
  HP2106201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2107901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2108901@meta.data[["Chemistry"]] <- "10Xv3"
  HP2110001@meta.data[["Chemistry"]] <- "10Xv3"
  HP2121601@meta.data[["Chemistry"]] <- "10Xv3"
  HP2123201@meta.data[["Chemistry"]] <- "10Xv3"
  HP2132801@meta.data[["Chemistry"]] <- "10Xv3"
  HP2202101@meta.data[["Chemistry"]] <- "10Xv3"
  
  
  HP2022801@meta.data[["Diabetes Status"]] <- "ND"
  SAMN15877725@meta.data[["Diabetes Status"]] <- "ND"
  HP2024001@meta.data[["Diabetes Status"]] <- "ND"
  HP2031401@meta.data[["Diabetes Status"]] <- "ND"
  HP2105501@meta.data[["Diabetes Status"]] <- "ND"
  HP2106201@meta.data[["Diabetes Status"]] <- "ND"
  HP2107001@meta.data[["Diabetes Status"]] <- "ND"
  HP2107901@meta.data[["Diabetes Status"]] <- "ND"
  HP2108601@meta.data[["Diabetes Status"]] <- "ND"
  HP2108901@meta.data[["Diabetes Status"]] <- "ND"
  HP2110001@meta.data[["Diabetes Status"]] <- "ND"
  HP2121601@meta.data[["Diabetes Status"]] <- "ND"
  HP2123201@meta.data[["Diabetes Status"]] <- "ND"
  HP2132801@meta.data[["Diabetes Status"]] <- "ND"
  HP2202101@meta.data[["Diabetes Status"]] <- "ND"
}

# Make list
pancreas_qadir <- list("HP2022801" = HP2022801, "SAMN15877725" = SAMN15877725, "HP2107001" = HP2107001, "HP2107901" = HP2107901,
                       "HP2024001" = HP2024001, "HP2105501" = HP2105501, "HP2108601" = HP2108601, "HP2108901" = HP2108901, 
                       "HP2031401" = HP2031401, "HP2110001" = HP2110001, "HP2123201" = HP2123201,
                       "HP2106201" = HP2106201, "HP2121601" = HP2121601, "HP2132801" = HP2132801, "HP2202101" = HP2202101
)

# Subset out all single cells
pancreas_qadir <- lapply(X = pancreas_qadir, 
                         FUN = function(x){
                           
                           # Subset out all singlets
                           Idents(x) <- "doublets"
                           x <- subset(x = x, idents = c("Singlet"))
                         })

#Save point
#saveRDS(pancreas_qadir, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_qadir.rds)")


# Load Processed data for analysis
#pancreas_subset <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_subset.rds)")
#pancreas_qadir <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_qadir.rds)")

# Join datasets
pancreas_combined <- c(pancreas_subset, pancreas_qadir)

#Save point
#saveRDS(pancreas_combined, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_combined.rds)")

# Load Processed data for analysis
#pancreas_combined <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_combined.rds)")

# merge data sets
sampleset <- c("HPAP-022", "HPAP_026", "HPAP_035", "HPAP_036", "HPAP_037", "HPAP_040", "HPAP_044", "HPAP_051", 
               "HPAP_052", "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", 
               "HPAP_063", "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", 
               "HPAP_081", "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099",
               "HPAP_100", "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109", "HP2022801", 
               "SAMN15877725", "HP2107001", "HP2107901", "HP2024001", "HP2105501", "HP2108601", "HP2108901", 
               "HP2031401", "HP2110001", "HP2123201", "HP2106201", "HP2121601", "HP2132801", "HP2202101")

pancreas_rna <- merge(pancreas_combined[[sampleset[[1]]]], 
                     y=pancreas_combined[sampleset[2:length(sampleset)]], 
                     project='pancreas', 
                     merge.data = TRUE)


# Inspect data
Idents(pancreas_rna) <- "Tissue Source"
VlnPlot(object = pancreas_rna, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3, pt.size = 0)

# Subset
pancreas_rna <- subset(x = pancreas_rna, subset = nFeature_RNA > 500 & nFeature_RNA < 8000 & percent.mt < 15)

#Save point
#saveRDS(pancreas_rna, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_rna.rds)")


# Load data
#pancreas_rna <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_rna.rds)")
#pancreas_combined <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_combined.rds)")


# Inspect data
pancreas_rna
table(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Sex"]])
table(pancreas_rna@meta.data[["Library"]])
table(pancreas_rna@meta.data[["Tissue Source"]])
table(pancreas_rna@meta.data[["Diabetes Status"]])
table(pancreas_rna@meta.data[["ancestry"]])
table(pancreas_rna@meta.data[["doublets"]])
table(pancreas_rna@meta.data[["Cell Type"]])

# #Split data on basis of disease status and calculate integration features on object list
# DefaultAssay(pancreas_rna) <- "SCT"
# integrationfeatures <- SelectIntegrationFeatures(pancreas_combined, nfeatures = 3000, verbose = TRUE)
# 
# # replacing variable features with integrationfeatures
# VariableFeatures(pancreas_rna, assay = "SCT") <- integrationfeatures
# VariableFeatures(pancreas_rna, assay = "RNA") <- integrationfeatures


#Perform basic threshold filtering and log normalization
DefaultAssay(pancreas_rna) <- "RNA"
pancreas_rna <- NormalizeData(pancreas_rna, normalization.method = "LogNormalize", scale.factor = 10000)
pancreas_rna <- FindVariableFeatures(pancreas_rna, selection.method = "vst", nfeatures = 2000)
pancreas_rna <- ScaleData(pancreas_rna, verbose = FALSE) %>% 
  RunPCA(pc.genes = pancreas_rna@assays$RNA@var.features, npcs = 20, verbose = FALSE)

# Re-level object@meta.data this just orders the actual metadata slot, so when you pull its already ordered
Idents(pancreas_rna) <- "ancestry"
pancreas_rna$ancestry_sex <- paste(Idents(pancreas_rna), pancreas_rna$Sex, sep = "_")
table(pancreas_rna@meta.data[["ancestry_sex"]])

#Run Harmony batch correction with library and tissue source covariates
Idents(pancreas_rna) <- "Library"
pancreas_rna <- RunHarmony(pancreas_rna, 
                           assay.use = "RNA",
                           reduction = "pca",
                           dims.use = 1:20,
                           group.by.vars = c('Library','Tissue Source', 'Chemistry'),
                           kmeans_init_nstart=20, kmeans_init_iter_max=100,
                           plot_convergence = TRUE)

unique(pancreas_rna@meta.data[["Library"]])
unique(pancreas_rna@meta.data[["Tissue Source"]])
unique(pancreas_rna@meta.data[["ancestry_sex"]])
unique(pancreas_rna@meta.data[["Cell Type"]])
unique(pancreas_rna@meta.data[["Chemistry"]])
table(pancreas_rna@meta.data[["Chemistry"]])

# Run UMAP
pancreas_rna <- RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)

# Clustering
pancreas_rna <- pancreas_rna %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 2.6, 2.8, 3, 3.2, 3.4, 3.6, 4, 4.6, 5, 6, 7, 8, 9, 10), method = 'igraph')

# Save file
saveRDS(pancreas_rna, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_rna.rds)")
})

# Cluster-tree analysis, looking appropriate non-anomalous clustering resolution
clustree(pancreas_rna, prefix = "RNA_snn_res.")

pancreas_rna <- FindClusters(pancreas_rna, algorithm=4,resolution = c(10), method = 'igraph')

# View clustering
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.4', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = pancreas_rna,
            features = c("SOX10"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

# Subclustering
Idents(pancreas_rna) <- "RNA_snn_res.6"
subset_clust <- subset(pancreas_rna, idents = c("1", "2", "3", "4", "5", "6", "7", "8", "9", "10",
                                                "11", "12", "13", "14", "15", "16", "17", "18", "19", "20",
                                                "21", "22", "23", "24", "25", "26", "27", "28", "29", "30", 
                                                "31", "32", "33", "34", "35", "36", "37", "38", "39", "40", 
                                                "41", "42", "43", "44", "45", "46", "47", "48", "49", "50", 
                                                "51", "52", "53", "54", "55", #"56", 
                                                "57", "58", "59", "60",
                                                "61", "62", "63", "64", "65", "66", #"67", 
                                                "68", "69", #"70", 
                                                "71", "72", "73", "74", "75", "76", "77", "78", "79", "80", 
                                                "81", "82", #"83", 
                                                "84", "85", "86", "87", "88", "89", "90" 
                                                ))

# Checking cluster loss
DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
DimPlot(subset_clust_test, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)

# As cells were subsetted re-run batch correction and cluster assignmnet
subset_clust <- RunHarmony(subset_clust, 
                           assay.use = "RNA",
                           reduction = "pca",
                           dims.use = 1:20,
                           group.by.vars = c('Library','Tissue Source', 'Chemistry'),
                           kmeans_init_nstart=20, kmeans_init_iter_max=100,
                           plot_convergence = TRUE)


# UMAP 
subset_clust <- RunUMAP(subset_clust, reduction = "harmony", dims = 1:20, return.model = TRUE)

#Neighbours + Clustering
subset_clust <- subset_clust %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1, 6), method = 'igraph')


DimPlot(subset_clust, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = subset_clust,
            features = c("HSPB1", "DNAJB6", "HSPH1", "GADD45B"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 100,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=TRUE)

DimPlot(pancreas_rna, reduction = 'umap', group.by = 'RNA_snn_res.6', label = TRUE, pt.size = 0.01, raster=FALSE)
FeaturePlot(object = pancreas_rna,
            features = c("MT-CO3"
            ),
            pt.size = 0.01,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=FALSE)

# Cluster assignment
table(subset_clust@meta.data$RNA_snn_res.6)
Idents(subset_clust) <- "RNA_snn_res.6"
subset_clust <- RenameIdents(subset_clust, 
                             "1" = "Alpha",
                             "2" = "Alpha",
                             "3" = "Alpha",
                             "4" = "Alpha",
                             "5" = "Alpha",
                             "6" = "Alpha",
                             "7" = "Acinar",
                             "8" = "Quiescent_Stellate",
                             "9" = "Alpha",
                             "10" = "Alpha+Beta",
                             "11" = "Acinar",
                             "12" = "Beta",
                             "13" = "Alpha",
                             "14" = "Alpha+Beta",
                             "15" = "Beta",
                             "16" = "Alpha",
                             "17" = "Ductal",
                             "18" = "Alpha",
                             "19" = "Delta",
                             "20" = "Ductal",
                             "21" = "Acinar",
                             "22" = "Alpha",
                             "23" = "Beta",
                             "24" = "Beta",
                             "25" = "Activated_Stellate",
                             "26" = "Ductal",
                             "27" = "Ductal",
                             "28" = "Acinar",
                             "29" = "Beta",
                             "30" = "Acinar",
                             "31" = "Alpha",
                             "32" = "Beta",
                             "33" = "Activated_Stellate",
                             "34" = "Alpha",
                             "35" = "Beta",
                             "36" = "Alpha",
                             "37" = "Endothelial",
                             "38" = "Beta",
                             "39" = "Beta",
                             "40" = "Endothelial",
                             "41" = "Acinar",
                             "42" = "Acinar",
                             "43" = "Activated_Stellate",
                             "44" = "Delta",
                             "45" = "Beta",
                             "46" = "Alpha",
                             "47" = "Beta",
                             "48" = "Alpha",
                             "49" = "Alpha+Beta",
                             "50" = "Alpha",
                             "51" = "Beta",
                             "52" = "Mast",
                             "53" = "Acinar",
                             "54" = "Alpha+Beta",
                             "55" = "Macrophage",
                             "56" = "Acinar",
                             "57" = "Acinar",
                             "58" = "Activated_Stellate",
                             "59" = "Beta",
                             "60" = "PPY",
                             "61" = "Acinar",
                             "62" = "PPY",
                             "63" = "Delta",
                             "64" = "Ductal",
                             "65" = "Endothelial",
                             "66" = "Quiescent_Stellate",
                             "67" = "Ductal",
                             "68" = "Acinar",
                             "69" = "Cycling_Endo",
                             "70" = "Alpha",
                             "71" = "Alpha+Beta",
                             "72" = "Beta",
                             "73" = "Endothelial",
                             "74" = "Acinar",
                             "75" = "Endothelial",
                             "76" = "Ductal",
                             "77" = "Acinar",
                             "78" = "Macrophage",
                             "79" = "Acinar",
                             "80" = "Schwann",
                             "81" = "Acinar",
                             "82" = "EndMT",
                             "83" = "Alpha",
                             "84" = "Beta",
                             "85" = "Macrophage",
                             "86" = "Acinar",
                             "87" = "Acinar"
                             )

# Check renaming
table(subset_clust@active.ident)

#plot <- DimPlot(subset_clust, reduction = "umap")
DefaultAssay(object = subset_clust) <- "RNA"
Idents(subset_clust, WhichCells(object = subset_clust, expression = GHRL > 50, slot = 'counts')) <- 'Epsilon'
DimPlot(subset_clust, reduction = "umap", label = TRUE)

# Subsetting Lymphocytes from Mast cells
# First subset all cells, they are called "Mast" in the primary object
lymphocytes <- subset(subset_clust, idents = "Mast")
table(lymphocytes@active.ident)
lymphocytes <- lymphocytes %>% 
  FindNeighbors(reduction = 'harmony', dims = 1:20) %>% 
  FindClusters(algorithm=4,resolution = c(0.5), method = 'igraph')

# Run UMAP
lymphocytes <- RunUMAP(lymphocytes, reduction = "harmony", dims = 1:20, return.model = TRUE)
DimPlot(lymphocytes, reduction = 'umap', label = FALSE, pt.size = 1, raster=FALSE)

# Cluster assignment
table(lymphocytes@meta.data$RNA_snn_res.0.5)
Idents(lymphocytes) <- "RNA_snn_res.0.5"
DimPlot(lymphocytes, reduction = 'umap', group.by = 'RNA_snn_res.0.5', label = TRUE, pt.size = 1, raster=FALSE)
FeaturePlot(object = lymphocytes,
            features = c("TPSAB1"
            ),
            pt.size = 1,
            cols = c("darkgrey", "red"),
            min.cutoff = 0,
            max.cutoff = 500,
            slot = 'counts',
            order = TRUE,
            raster=TRUE)
table(lymphocytes@active.ident)
lymphocytes <- RenameIdents(lymphocytes, 
                            "1" = "Lymphocyte",
                            "2" = "Mast",
                            "3" = "Mast",
                            "4" = "Mast"
)
table(lymphocytes@active.ident)

# Generate a new column called celltype_qadir in the metadata copying all Ident info there, this is active idents so check
table(subset_clust@active.ident)
subset_clust$celltype_qadir <- as.character(Idents(subset_clust))
table(subset_clust@meta.data$celltype_qadir)

# Change the information of cells containing sub-cluster information
subset_clust$celltype_qadir[Cells(lymphocytes)] <- paste("c1",Idents(lymphocytes))
DimPlot(pbmc_small, group.by = "sub_cluster")


# Saving this information in the metadata slot
table(Idents(pancreas.combined.h.s))
pancreas.combined.h.s$celltype <- Idents(pancreas.combined.h.s)
head(pancreas.combined.h.s@meta.data)
# 
#
#
#
#
#
#
#
#
#
#
#
#
# View clustering
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
pancreas_rna <- FindClusters(object = pancreas_rna, algorithm=4, resolution = 0, method = 'igraph')


Idents(pancreas_rna) <- "Chemistry"

# QC for Ruth split on chemistry
v2 <- subset(pancreas_rna, idents = c("10Xv2"))
v3 <- subset(pancreas_rna, idents = c("10Xv3"))

Idents(pancreas_rna) <- "Chemistry"
DimPlot(pancreas_rna, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)

p1 <- DimPlot(v2, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
p2 <- DimPlot(v3, reduction = 'umap', label = FALSE, pt.size = 0.01, raster=FALSE)
p1|p2

RunUMAP(pancreas_rna, reduction = "harmony", dims = 1:20, return.model = TRUE, reduction.name = 'umap')
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)
DimPlot(pancreas_rna, reduction = 'harmony', label = FALSE, pt.size = 0.01, raster=FALSE)









