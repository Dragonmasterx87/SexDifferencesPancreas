# NAalysis of Ruth's data
# Hpap dataset
hpap <- readRDS(file=r"(C:\Users\mqadir\Downloads\hpap_islet_scRNAseq.rds)")
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

pancreas.list <- lapply(X = pancreas.list, 
                        FUN = function(x){
                          x[['percent.mt']] <- PercentageFeatureSet(x, 
                                                                    pattern = '^MT-')
                          x <- SCTransform(x,
                                           vars.to.regress = "percent.mt", 
                                           verbose = TRUE, 
                                           return.only.var.genes = FALSE, 
                                           vst.flavor = "v2")
                          })

# Save file because it takes ages to run
# saveRDS(pancreas.list, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_list.rds)")
# pancreas.list <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_list.rds)")

# Pulling data for analysis
pancreas.list.use <- pancreas.list[c("HPAP_019", "HPAP-022", "HPAP_024", "HPAP_026", "HPAP_029", "HPAP_035", "HPAP_036", "HPAP_037",
                                     "HPAP_040", "HPAP_043", "HPAP_044", "HPAP_045", "HPAP_049", "HPAP_050", "HPAP_051", "HPAP_052",
                                     "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", "HPAP_063",
                                     "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", "HPAP_081",
                                     "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099", "HPAP_100", 
                                     "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109")
                                   ]

#Save point
#saveRDS(pancreas.list.use, file = r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_list_subset.rds)")

{
  HP2022801 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2022801.rds)")
  SAMN15877725 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\SAMN15877725.rds)")
  HP2024001 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2024001.rds)")
  HP2031401 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2031401.rds)")
  HP2105501 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2105501.rds)")
  HP2106201 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2106201.rds)")
  HP2107001 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2107001.rds)")
  HP2107901 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2107901.rds)")
  HP2108601 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2108601.rds)")
  HP2108901 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2108901.rds)")
  HP2110001 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2110001.rds)")
  HP2121601 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2121601.rds)")
  HP2123201 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2123201.rds)")
  HP2132801 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2132801.rds)")
  HP2202101 <- readRDS(r"(D:\3.CodingScriptsandData\SexBasedstudy\RDS files\rna\CoreArch_CurrentR\HP Files after doubletfinder\HP2202101.rds)")
}

# Adjust and add metadata
# Sex specific Metadata addition
{
  HP2022801$sex <- "M"
  SAMN15877725$sex <- "M"
  HP2024001$sex <- "F"
  HP2031401$sex <- "M"
  HP2105501$sex <- "F"
  HP2106201$sex <- "F"
  HP2107001$sex <- "M"
  HP2107901$sex <- "M"
  HP2108601$sex <- "F"
  HP2108901$sex <- "F"
  HP2110001$sex <- "M"
  HP2121601$sex <- "F"
  HP2123201$sex <- "M"
  HP2132801$sex <- "F"
  HP2202101$sex <- "F"


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
pancreas_subset <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_list_subset.rds)")
pancreas_qadir <- readRDS(r"(E:\2.SexbasedStudyCurrent\RDS files\Ruth_data\pancreas.list\pancreas_qadir.rds)")

# Join datasets
pancreas_combined <- c(pancreas_subset, pancreas_qadir)

# merge data sets
sampleset <- c("HPAP_019", "HPAP-022", "HPAP_024", "HPAP_026", "HPAP_029", "HPAP_035", "HPAP_036", "HPAP_037",
               "HPAP_040", "HPAP_043", "HPAP_044", "HPAP_045", "HPAP_049", "HPAP_050", "HPAP_051", "HPAP_052",
               "HPAP_053", "HPAP_054", "HPAP_056", "HPAP_057", "HPAP_058", "HPAP_059", "HPAP_061", "HPAP_063",
               "HPAP_065", "HPAP_070", "HPAP_074", "HPAP_075", "HPAP_077", "HPAP_079", "HPAP_080", "HPAP_081",
               "HPAP_082", "HPAP_083", "HPAP_085", "HPAP_088", "HPAP_090", "HPAP_091", "HPAP_099", "HPAP_100", 
               "HPAP_101", "HPAP_103", "HPAP_105", "HPAP_106", "HPAP_108", "HPAP_109", "HP2022801", "SAMN15877725",
               "HP2107001", "HP2107901", "HP2024001", "HP2105501", "HP2108601", "HP2108901", "HP2031401", "HP2110001",
               "HP2123201", "HP2106201", "HP2121601", "HP2132801", "HP2202101")

merged_data <- merge(pancreas_combined[[sampleset[[1]]]], 
                     y=pancreas_combined[sampleset[2:length(sampleset)]], 
                     project='pancreas', 
                     merge.data = TRUE)

# Inspect data
merged_data
Idents(merged_data) <- "Tissue Source"
VlnPlot(object = pancreas.combined.h, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), ncol = 3)

















