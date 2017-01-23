BINDIR := bin
TARGETS := $(addprefix $(BINDIR)/,abawaca abawaca-build)
_dummy := $(shell mkdir -p $(BINDIR))

.PHONY: all

all: $(TARGETS)

$(BINDIR)/abawaca: $(BINDIR)/abawaca.o $(BINDIR)/tinymt32.o $(BINDIR)/ScafDpData.o $(BINDIR)/Dimension.o $(BINDIR)/ClusterData.o $(BINDIR)/ClusterWriter.o $(BINDIR)/Cluster.o $(BINDIR)/misc.o $(BINDIR)/ClusterSeparator.o $(BINDIR)/ClusterSeparatorBySensitivitySpecificity.o $(BINDIR)/ClusterSeparatorSplitScafs.o $(BINDIR)/SCGdb.o $(BINDIR)/String.o $(BINDIR)/common.o $(BINDIR)/ClusterQuality.o
	g++ -o $(BINDIR)/abawaca $(BINDIR)/abawaca.o $(BINDIR)/tinymt32.o $(BINDIR)/ScafDpData.o $(BINDIR)/Dimension.o $(BINDIR)/ClusterData.o $(BINDIR)/ClusterWriter.o $(BINDIR)/Cluster.o $(BINDIR)/misc.o $(BINDIR)/ClusterSeparator.o $(BINDIR)/ClusterSeparatorBySensitivitySpecificity.o $(BINDIR)/ClusterSeparatorSplitScafs.o $(BINDIR)/SCGdb.o $(BINDIR)/String.o $(BINDIR)/common.o $(BINDIR)/ClusterQuality.o -pthread -std=c++11

$(BINDIR)/abawaca-build: $(BINDIR)/ReadMapping.o $(BINDIR)/ReadMappingReader.o $(BINDIR)/abawaca-build.o $(BINDIR)/common.o $(BINDIR)/String.o
	g++ -o $(BINDIR)/abawaca-build $(BINDIR)/ReadMapping.o $(BINDIR)/ReadMappingReader.o $(BINDIR)/abawaca-build.o $(BINDIR)/common.o $(BINDIR)/String.o -pthread -std=c++11

$(BINDIR)/abawaca-build.o: src/abawaca-build.cpp src/ReadMappingReader.h src/common.h src/ReadMapping.h src/String.h src/Sequence.h src/bio_exceptions.h src/SeqIO.h src/SeqIORead_fasta.h src/SeqIORead.h
	g++ -c src/abawaca-build.cpp -o $(BINDIR)/abawaca-build.o -pthread -std=c++11

$(BINDIR)/abawaca.o: src/abawaca.cpp src/misc.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/ClusterData.h src/Dimension.h src/ClusterSeparatorBySensitivitySpecificity.h src/ClusterSeparator.h src/Cluster.h src/Semaphore.h src/SCGdb.h src/ClusterQuality.h src/ClusterSeparatorSplitScafs.h src/ClusterWriter.h src/SeqIOWrite_fasta.h src/SeqIOWrite.h src/SeqIO.h
	g++ -c src/abawaca.cpp -o $(BINDIR)/abawaca.o -pthread -std=c++11

$(BINDIR)/Cluster.o: src/Cluster.cpp src/Cluster.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/Dimension.h
	g++ -c src/Cluster.cpp -o $(BINDIR)/Cluster.o -pthread -std=c++11

$(BINDIR)/ClusterData.o: src/ClusterData.cpp src/ClusterData.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/Dimension.h
	g++ -c src/ClusterData.cpp -o $(BINDIR)/ClusterData.o -pthread -std=c++11

$(BINDIR)/ClusterQuality.o: src/ClusterQuality.cpp src/ClusterQuality.h src/Cluster.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/Dimension.h src/SCGdb.h
	g++ -c src/ClusterQuality.cpp -o $(BINDIR)/ClusterQuality.o -pthread -std=c++11

$(BINDIR)/ClusterSeparatorBySensitivitySpecificity.o: src/ClusterSeparatorBySensitivitySpecificity.cpp src/ClusterSeparatorBySensitivitySpecificity.h src/ClusterSeparator.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/ClusterData.h src/Dimension.h src/Cluster.h src/Semaphore.h src/SCGdb.h src/ClusterQuality.h
	g++ -c src/ClusterSeparatorBySensitivitySpecificity.cpp -o $(BINDIR)/ClusterSeparatorBySensitivitySpecificity.o -pthread -std=c++11

$(BINDIR)/ClusterSeparator.o: src/ClusterSeparator.cpp src/ClusterSeparator.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/ClusterData.h src/Dimension.h src/Cluster.h src/Semaphore.h src/SCGdb.h src/ClusterQuality.h src/misc.h
	g++ -c src/ClusterSeparator.cpp -o $(BINDIR)/ClusterSeparator.o -pthread -std=c++11

$(BINDIR)/ClusterSeparatorSplitScafs.o: src/ClusterSeparatorSplitScafs.cpp src/ClusterSeparatorSplitScafs.h src/ClusterSeparator.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/ClusterData.h src/Dimension.h src/Cluster.h src/Semaphore.h src/SCGdb.h src/ClusterQuality.h
	g++ -c src/ClusterSeparatorSplitScafs.cpp -o $(BINDIR)/ClusterSeparatorSplitScafs.o -pthread -std=c++11

$(BINDIR)/ClusterWriter.o: src/ClusterWriter.cpp src/ClusterWriter.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/ClusterData.h src/Dimension.h
	g++ -c src/ClusterWriter.cpp -o $(BINDIR)/ClusterWriter.o -pthread -std=c++11

$(BINDIR)/common.o: src/common.cpp src/common.h
	g++ -c src/common.cpp -o $(BINDIR)/common.o -pthread -std=c++11

$(BINDIR)/Dimension.o: src/Dimension.cpp src/Dimension.h
	g++ -c src/Dimension.cpp -o $(BINDIR)/Dimension.o -pthread -std=c++11

$(BINDIR)/misc.o: src/misc.cpp src/misc.h src/TinyMT-src-1.0.1/tinymt/tinymt32.h
	g++ -c src/misc.cpp -o $(BINDIR)/misc.o -pthread -std=c++11

$(BINDIR)/ReadMapping.o: src/ReadMapping.cpp src/ReadMapping.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h
	g++ -c src/ReadMapping.cpp -o $(BINDIR)/ReadMapping.o -pthread -std=c++11

$(BINDIR)/ReadMappingReader.o: src/ReadMappingReader.cpp src/ReadMappingReader.h src/common.h src/ReadMapping.h src/String.h src/Sequence.h src/bio_exceptions.h src/SeqIO.h
	g++ -c src/ReadMappingReader.cpp -o $(BINDIR)/ReadMappingReader.o -pthread -std=c++11

$(BINDIR)/ScafDpData.o: src/ScafDpData.cpp src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h src/SeqIORead_fasta.h src/SeqIORead.h src/SeqIO.h
	g++ -c src/ScafDpData.cpp -o $(BINDIR)/ScafDpData.o -pthread -std=c++11

$(BINDIR)/SCGdb.o: src/SCGdb.cpp src/SCGdb.h src/ScafDpData.h src/common.h src/String.h src/Sequence.h src/bio_exceptions.h
	g++ -c src/SCGdb.cpp -o $(BINDIR)/SCGdb.o -pthread -std=c++11

$(BINDIR)/String.o: src/String.cpp src/String.h src/common.h src/Sequence.h src/bio_exceptions.h
	g++ -c src/String.cpp -o $(BINDIR)/String.o -pthread -std=c++11

$(BINDIR)/tinymt32.o: src/TinyMT-src-1.0.1/tinymt/tinymt32.c src/TinyMT-src-1.0.1/tinymt/tinymt32.h
	g++ -c src/TinyMT-src-1.0.1/tinymt/tinymt32.c -o $(BINDIR)/tinymt32.o -pthread -std=c++11

$(BINDIR):
	mkdir -p $(BINDIR)
