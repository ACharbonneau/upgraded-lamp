# There are no unique identifiers in this project until the DNA prep,
# and even in the DNA prep, the unique identifiers aren't unique because
# ~80 samples were sequenced twice.
# One of these samples, 12059, was NOT supposed to be sequenced twice, but appears
# twice in the sequencing data, in the same plate. Since I have the correct number
# of samples, this means another sample must not have been sequenced. Or at least
# was mislabeled. This sample is 13245, which was supposed to be sequenced twice,
# but wasn't. I asked Jeff to check his metadata collection to see if he could
# determine whether 12059 really was sequenced twice, or if one of them is 13245
# just mislabeled.

# Jeff checked and now thinks that:
# "the sample in Cell C4 of plate 4 was really 12054, not 12059, so it is KXF2(2)/102."

# Since the two samples labeled 12059 are
# in the same plate, they get the same unique ID, so I've added a really stupid bit
# of code to add the sequencing well to the unique ID as well.

# Furthermore the ID names change subtley between datasheets.
# This script just munges the metadata into one useable format with unique IDs.


rm(list = ls())

# Install function for packages
packages<-function(x){
  x<-as.character(match.call()[[2]])
  if (!require(x,character.only=TRUE)){
    install.packages(pkgs=x,repos="http://cran.r-project.org")
    require(x,character.only=TRUE)
  }
}
packages(tidyr)
packages(plyr)
packages(dplyr)


## Get global data frames

DNA_data <- read.csv("../Metadata/OriginalFiles/MetadataAll.txt", sep = "\t", header = F,
                     colClasses = "factor")

colnames(DNA_data) <- c("ID", "Type_Year", "Species", "Cross", "Indiv", "Date", "Prep")

# File from Jeff has all of these wrong :(
DNA_data$Species[DNA_data$Cross=="SAES"] <- "Rrl"
DNA_data$Species[DNA_data$Cross=="RA808"] <- "Rrl"
DNA_data$Species[DNA_data$Cross=="YEIL_CLNC"] <- "Rros"
DNA_data$Cross[DNA_data$Cross=="SPEU"] <- "SPNK"

DNA_data$Cross <- revalue(DNA_data$Cross, c("YEIL_CLNC" = "YEIL"))

DNA_data <- droplevels(DNA_data)

Pedigree <- read.csv("../Metadata/OriginalFiles/Pedigree.csv", colClasses = "factor")


## They sequenced lots of samples twice, but used the same names, which breaks
## the deconvolution step. This loop gives each sample a unique name and writes
## out a new sequence list file for each lane with the unique names

Geno_dir <- file.path("../Metadata/PlateInfoSeq/")

Geno_list <- list.files(Geno_dir, pattern = "*gz.keys.txt")

All_geno_data <- list()

for( X in 1:length(Geno_list)){
  TempFile <- read.table(paste(Geno_dir, Geno_list[ X ], sep=""), sep = "\t", header = T, na.strings = "")
  TempFile$UniqID <- paste(TempFile$DNASample, gsub( '.*_([0-9]{1,2})', "\\1", TempFile$LibraryPlate, ignore.case = FALSE ), TempFile$SampleDNA_Well, "F", sep=".")
  ComboTemp <- left_join(TempFile, DNA_data, by=c("DNASample"="ID"))
  All_geno_data[[X]] <- ComboTemp
  write.table(select(ComboTemp, -Pedigree, -Population, -SeedLot), paste(Geno_dir, substr(Geno_list[X], 1, 11), ".unique.txt", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
  write.table(cbind(as.character(ComboTemp$Barcode), paste(ComboTemp$Cross, ComboTemp$UniqID, sep="_")), paste(Geno_dir, substr(Geno_list[X], 1, 11), "_fastq.gz.barcodes", sep=""), sep="\t", row.names=F, col.names=F, quote=F)
}

rm(TempFile, ComboTemp)

All_geno_data <- tbl_df( do.call( "rbind", All_geno_data ))



# Get F0 data and combine

F0_Pheno_data <- read.csv("../Metadata/OriginalFiles/QTLParentalMeasurements.csv", colClasses =
                            c(rep("factor", 4), rep("numeric", 6), "factor"))

F0_Pheno_data$Indiv <- paste(F0_Pheno_data$Matriline, "/", F0_Pheno_data$Offspr, sep="")

AE_F0_DNA <- filter( DNA_data, Type_Year=="AS_EX_QTL_Parents_2011")

AE_F0_DNA <- droplevels(AE_F0_DNA)

AE_F0_DNA <- left_join(AE_F0_DNA, F0_Pheno_data)

write.csv(x = AE_F0_DNA, file= "../Metadata/AE_F0_merge.csv")

# Get F1 data and combine

F1_Pheno_data <- read.csv("../Metadata/OriginalFiles/QTLF1measurements.csv", colClasses =
                            c(rep("factor", 4), rep("numeric", 6), "factor"))

F1_Pheno_data$UID <- paste(F1_Pheno_data$Plant.ID, "/", F1_Pheno_data$Offspr, sep="")

AE_F1_DNA <- filter( DNA_data, Type_Year=="AS_EX_QTL_F1s_2011")

AE_F1_DNA$UID <- paste(AE_F1_DNA$Cross, AE_F1_DNA$Indiv, sep="")

AE_F1_DNA <- left_join(AE_F1_DNA, F1_Pheno_data)

AE_F1_DNA <- droplevels(AE_F1_DNA)

write.csv(x = AE_F1_DNA, file = "../Metadata/AE_F1_merge.csv")


# Get F2 data and combine

F2_Pheno_data <- read.csv("../Metadata/OriginalFiles/Exsertion_F2_F2s.csv", colClasses =
                         c(rep("factor", 7), rep("numeric", 7), "factor"))

colnames(F2_Pheno_data)[1] <- "CrossX"

## Filter out F2s from Parents from F1s

AE_F2_DNA <- filter( DNA_data, Type_Year != "SigSelection" &
                       Type_Year!="AS_EX_QTL_Parents_2011" &
                       Type_Year!="AS_EX_QTL_F1s_2011") %>%
              separate( Indiv, into = c("Family", "F2" ), sep = "/")

AE_F2_DNA <- droplevels(AE_F2_DNA)

AE_F2_DNA$F2 <- as.factor(AE_F2_DNA$F2)

## Merge F2s sequencing and phenotype data

CrossList <- c(KF2="KxK", KXF2="KxR", RF2="RxR", RXF2="RxK")

AE_F2_DNA$Family <- factor(gsub( "[()]", "", AE_F2_DNA$Family ))

AE_F2_DNA$CrossX <- as.factor(CrossList[AE_F2_DNA$Cross])

AE_F2_DNA <- left_join(AE_F2_DNA, F2_Pheno_data)

AE_F2_DNA$F2s.produced <- as.factor(paste(AE_F2_DNA$Cross, "(", AE_F2_DNA$Family, ")", sep = ""))

AE_F2_DNA <- left_join(AE_F2_DNA, Pedigree)

write.csv(x = AE_F2_DNA, file = "../Metadata/AE_F2_merge.csv")


# Write out metadata for running subsets
# I need a file for each subset that I can feed into samtools or similar:
# a file of filenames, one file per line


# This is the filenames with two columns of metatdata

ForStacksAE <- rbind(select(AE_F2_DNA, ID, Cross, Type_Year),
      select(AE_F0_DNA, ID, Cross, Type_Year),
      select(AE_F1_DNA, ID, Cross, Type_Year))

ForStacksAEUniq <- dplyr::left_join(ForStacksAE, All_geno_data, by=c("ID"="DNASample", "Cross"="Cross", "Type_Year"="Type_Year"))

ForStacksAEUniq <- unique(ForStacksAEUniq)

ForStacksAEUniq$UniqID <- paste(ForStacksAEUniq$UniqID, ".trimmed.fq_q20.final.bam", sep = "")

write.table(x = select(ForStacksAEUniq, UniqID, Cross, Type_Year), file = "../Metadata/AE_data.pop",
            quote = F, sep = "\t", col.names = F, row.names = F)


Just_F0s <- ForStacksAEUniq$UniqID[ForStacksAEUniq$Cross == "KH" |
                                     ForStacksAEUniq$Cross == "KL" |
                                     ForStacksAEUniq$Cross == "RH" |
                                     ForStacksAEUniq$Cross == "RL" ]

Just_F1s <- ForStacksAEUniq$UniqID[ForStacksAEUniq$Cross == "KF1A" |
                                     ForStacksAEUniq$Cross == "KF1B" |
                                     ForStacksAEUniq$Cross == "KXF1A" |
                                     ForStacksAEUniq$Cross == "KXF1B" |
                                     ForStacksAEUniq$Cross == "RF1A" |
                                     ForStacksAEUniq$Cross == "RF1B" |
                                     ForStacksAEUniq$Cross == "RXF1A" |
                                     ForStacksAEUniq$Cross == "RXF1B" ]

Just_F2s <- ForStacksAEUniq$UniqID[ForStacksAEUniq$Cross == "KF2" |
                                     ForStacksAEUniq$Cross == "KXF2" |
                                     ForStacksAEUniq$Cross == "RF2" |
                                     ForStacksAEUniq$Cross == "RXF2" ]


#These would make a single line of individuals with a flag between

#write.table(paste( " -s ./", Just_F0s, sep=""), file = "../Metadata/AE_F0_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
#write.table(paste( " -s ./", Just_F1s, sep=""), file = "../Metadata/AE_F1_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
#write.table(paste( " -s ./", Just_F2s, sep=""), file = "../Metadata/AE_F2_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")

#write.table(paste( " -s ./", Just_F0s, sep=""), file = "../Metadata/AE_Mapping_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
#write.table(paste( " -s ./", Just_F2s, sep=""), file = "../Metadata/AE_Mapping_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "", append = T)


#write.table(paste( " -p ./", Just_F0s, sep=""), file = "../Metadata/AE_Mapping_cs_stacks_geno", quote = F, col.names = F, row.names = F, eol = "")
#write.table(paste( " -r ./", Just_F2s, sep=""), file = "../Metadata/AE_Mapping_cs_stacks_geno", quote = F, col.names = F, row.names = F, eol = "", append = T)

#These are the files for samtools

write.table(Just_F0s, file = "../Metadata/AE_F0_list", quote = F, col.names = F, row.names = F)
write.table(Just_F1s, file = "../Metadata/AE_F1_list", quote = F, col.names = F, row.names = F)
write.table(Just_F2s, file = "../Metadata/AE_F2_list", quote = F, col.names = F, row.names = F)

write.table(Just_F0s, file = "../Metadata/AE_Mapping_list", quote = F, col.names = F, row.names = F)
write.table(Just_F2s, file = "../Metadata/AE_Mapping_list", quote = F, col.names = F, row.names = F, append = T)


## STACKS for signatures of selection

ForStacksSS <- filter( DNA_data, Type_Year == "SigSelection") %>%
  select(ID, Cross, Species)

ForStacksSSUniq <- dplyr::left_join(ForStacksSS, All_geno_data, by=c("ID"="DNASample", "Cross"="Cross"))
ForStacksSSUniq <- droplevels(ForStacksSSUniq)

ForStacksSSUniq$UniqID <- paste(ForStacksSSUniq$UniqID, ".trimmed.fq_q20.final.bam", sep = "")

write.table(x = select(ForStacksSSUniq, UniqID, Cross, Species.x), file = "../Metadata/SS_data.pop",
            quote = F, sep = "\t", col.names = F, row.names = F)

# Just RRR data
Just_Rrr <- ForStacksSSUniq$UniqID[ForStacksSSUniq$Species.x == "Rrr"]
#write.table(paste( " -s ./", Just_Rrr, sep=""), file = "../Metadata/SS_Rrr_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
write.table(Just_Rrr, file = "../Metadata/SS_Rrr_list", quote = F, col.names = F, row.names = F)

#Just landra
Just_landra <- ForStacksSSUniq$UniqID[ForStacksSSUniq$Species.x == "Rrl"]
#write.table(paste( " -s ./", Just_landra, sep=""), file = "../Metadata/SS_Rrl_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
write.table(Just_landra, file = "../Metadata/SS_Rrl_list", quote = F, col.names = F, row.names = F)

#Just crops
Just_Rsat <- ForStacksSSUniq$UniqID[ForStacksSSUniq$Species.x == "Rsat"]
#write.table(paste( " -s ./", Just_Rsat, sep=""), file = "../Metadata/SS_Rsat_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
write.table(Just_Rsat, file = "../Metadata/SS_Rsat_list", quote = F, col.names = F, row.names = F)

#RRR and landra
Rrr_landra <- ForStacksSSUniq$UniqID[ForStacksSSUniq$Species.x == "Rrr" | ForStacksSSUniq$Species.x == "Rrl"]
#write.table(paste( " -s ./", Rrr_landra, sep=""), file = "../Metadata/SS_RrrRrl_cs_stacks_list", quote = F, col.names = F, row.names = F, eol = "")
write.table(Rrr_landra, file = "../Metadata/SS_RrrRrl_list", quote = F, col.names = F, row.names = F)




#All signature of selection
All_SS <- ForStacksSSUniq$UniqID
write.table(All_SS, file = "../Metadata/SS_list", quote = F, col.names = F, row.names = F)




# Write out ChooseSigSel.sh

AEcommand <- cbind("cp", select(ForStacksAEUniq, UniqID), "AE_Deconvoluted/", select(ForStacksAEUniq, UniqID))
AEcommand$Folder <- paste(AEcommand[,3], AEcommand[,4], ".sam", sep = "")
colnames(AEcommand) <- c("cp", "UniqID", "JustFolder", "UniqRepeat", "Folder")


SScommand <- cbind("cp", select(ForStacksSSUniq, UniqID), "SigSelection/", select(ForStacksSSUniq, UniqID))
SScommand$Folder <- paste(SScommand[,3], SScommand[,4], ".sam", sep = "")
colnames(SScommand) <- c("cp", "UniqID", "JustFolder", "UniqRepeat", "Folder")
command <- rbind(AEcommand, SScommand)
command <- mutate(command, UniqID=paste("*",UniqID, ".fq.sam", sep=""))
command <- select(command, cp, UniqID, Folder)

#write.table(x = command, file = "1.4_ChooseSigSel.sh", quote = F,
#            sep = " ", col.names = F, row.names = F)

# Write out full data with phenotypes for SigSelection

FullSS <- filter( DNA_data, Type_Year == "SigSelection")
FullSSUniq <- dplyr::left_join(FullSS, All_geno_data, by=c("ID"="DNASample", "Cross"="Cross", "Type_Year"="Type_Year"))
FullSSUniq <- droplevels(FullSSUniq)

GeographyList <- c(AFFR="France", AROL="NA", BINY="nonNative", DEES="Spain", ESNK="NA",
                   GMIL="pugiformis", MAES="Spain", NAAU="nonNative", NELO="NA",
                   OIBG="NA", PBFR="France", RA226="Algeria", RA444="Italy",
                   RA761="Turkey", RA808="Turkey", RABG="NA", RACA="NA", SAES="Spain",
                   SPNK="NA", TOBG="NA", YEIL="pugiformis")

TaxonomyList <- c(AFFR="raphanistrum", AROL="oleifera", BINY="raphanistrum", DEES="raphanistrum",
                  ESNK="european", GMIL="pugiformis", MAES="raphanistrum", NAAU="raphanistrum", NELO="daikon",
                  OIBG="oleifera", PBFR="landra", RA226="landra", RA444="landra",
                  RA761="landra", RA808="landra", RABG="caudatus", RACA="caudatus", SAES="landra",
                  SPNK="european", TOBG="daikon", YEIL="pugiformis")


HabitatList <- c(AFFR="agricultural", AROL="cultivar", BINY="agricultural", DEES="natural",
                 ESNK="cultivar", GMIL="natural", MAES="disturbed", NAAU="agricultural", NELO="cultivar",
                 OIBG="cultivar", PBFR="natural", RA226="natural", RA444="natural",
                 RA761="natural", RA808="natural", RABG="cultivar", RACA="cultivar", SAES="natural",
                 SPNK="cultivar", TOBG="cultivar", YEIL="natural")

LocationList <- c(AFFR="raphNatW", AROL="oleifera", BINY="raphNN", DEES="raphNatW",
                  ESNK="european", GMIL="pugiformis", MAES="raphNatW", NAAU="raphNN", NELO="daikon",
                  OIBG="oleifera", PBFR="landra", RA226="landra", RA444="landra",
                  RA761="landra", RA808="landra", RABG="caudatus", RACA="caudatus", SAES="landra",
                  SPNK="european", TOBG="daikon", YEIL="pugiformis")

SpeciesList <- c(AFFR="raphanistrum", AROL="sativus", BINY="raphanistrum", DEES="raphanistrum",
                 ESNK="sativus", GMIL="pugiformis", MAES="raphanistrum", NAAU="raphanistrum", NELO="sativus",
                 OIBG="sativus", PBFR="raphanistrum", RA226="raphanistrum", RA444="raphanistrum",
                 RA761="raphanistrum", RA808="raphanistrum", RABG="sativus", RACA="sativus", SAES="raphanistrum",
                 SPNK="sativus", TOBG="sativus", YEIL="pugiformis")



OrderList <- c(AFFR=1, AROL=7, BINY=2, DEES=3,
               ESNK=1, GMIL=1, MAES=4, NAAU=5, NELO=3,
               OIBG=8, PBFR=1, RA226=3, RA444=4,
               RA761=5, RA808=6, RABG=5, RACA=6, SAES=2,
               SPNK=2, TOBG=4, YEIL=2)

FullSSUniq$Geo <- as.factor(GeographyList[FullSSUniq$Cross])
FullSSUniq$Taxon <- as.factor(TaxonomyList[FullSSUniq$Cross])
FullSSUniq$Habit <- as.factor(HabitatList[FullSSUniq$Cross])
FullSSUniq$locals <- as.factor(LocationList[FullSSUniq$Cross])
FullSSUniq$Species <- as.factor(SpeciesList[FullSSUniq$Cross])
FullSSUniq$Order <- as.factor(OrderList[FullSSUniq$Cross])
FullSSUniq$Order <- as.numeric(FullSSUniq$Order)

FullSSUniq <- select(FullSSUniq, -Pedigree, -Population, -SeedLot, -Species.y)

colnames(FullSSUniq) <- c("ID","Type_Year","STACKSspecies","Pop","Indiv","Date",
                               "Prep","Flowcell","Lane","Barcode","LibraryPlate",
                               "Row","Col","LibraryPrepID","LibraryPlateID","Enzyme",
                               "BarcodeWell","DNA_Plate","SampleDNA_Well","Genus",
                               "FullSampleName","UniqID","Geo","Taxon","Habit","locals",
                               "Species","Order")

write.csv(FullSSUniq, "../Metadata/SigSelectionMeta.csv")

