#!/usr/bin/env Rscript

set.seed(seed=19791123)

# Run AlphaSim
cat(19791123, file="Seed.txt")
system(command="./AlphaSim")

# Remove seed
system(command="rm -f Seed.txt SeedOld.txt")

# We do not need data on Selection decisions
system(command="rm -rf Selection")

# Keep only some relevant files
# ... pedigree
system(command="tail -n 400 SimulatedData/Pedigree.txt | awk '{ print $3\" \"$4\" \"$5 }' > ExamplePedigree.txt")
system(command="tail -n 400 SimulatedData/Pedigree.txt | awk '{ print $3\" \"$1        }' > ExampleGeneration.txt")
# ... marker array
system(command="mv -f SimulatedData/AllIndividualsSnpChips/Chip1Genotype.txt ExampleArrayGenotype.txt")
system(command="mv -f SimulatedData/AllIndividualsSnpChips/Chip1Phase.txt    ExampleArrayHaplotype.txt")
system(command="mv -f SimulatedData/Chip1SnpInformation.txt                  ExampleArrayMap.txt")
# ... sequence
system(command="mv -f SimulatedData/FullSequenceGenotype.txt            ExampleSequenceGenotype.txt")
system(command="mv -f SimulatedData/FullSequencePhase.txt               ExampleSequenceHaplotype.txt")
system(command="mv -f SimulatedData/FullSequencePhysicalMapIncluded.txt ExampleSequenceMap.txt")
system(command="mv -f SimulatedData/NbOfSegSitesPerChrom.txt            ExampleSequenceNSegSites.txt")
# ... cleanup
system(command="rm -rf SimulatedData")

# TODO: recombination tracking etc
system(command="rm -rf Chromosomes")

# A debug file?
system(command="rm -f temp")

# Make an incomplete pedigree
Ped <- read.table(file="ExamplePedigree.txt")
nInd <- nrow(Ped)
SelFather <- sample(1:nInd, size=nInd*0.10)
Ped[[2]][SelFather] <- 0
SelMother <- sample(1:nInd, size=nInd*0.20)
Ped[[3]][SelMother] <- 0
write.table(x=Ped, file="ExamplePedigreeIncomplete.txt", row.names=FALSE, col.names=FALSE)
