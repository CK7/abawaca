abawaca: main.o tinymt32.o ScafDpData.o Dimension.o ClusterData.o ClusterSeparator.o ClusterWriter.o Cluster.o misc.o
	g++ -o abawaca main.o tinymt32.o ScafDpData.o Dimension.o ClusterData.o ClusterSeparator.o ClusterWriter.o Cluster.o misc.o

ClusterSeparator.o: src/ClusterSeparator.cpp
	g++ -c src/ClusterSeparator.cpp -o ClusterSeparator.o

Cluster.o: src/Cluster.cpp
	g++ -c src/Cluster.cpp -o Cluster.o

ClusterWriter.o: src/ClusterWriter.cpp
	g++ -c src/ClusterWriter.cpp -o ClusterWriter.o

main.o: src/main.cpp
	g++ -c src/main.cpp -o main.o   

ClusterData.o: src/ClusterData.cpp
	g++ -c src/ClusterData.cpp -o ClusterData.o

Dimension.o: src/Dimension.cpp
	g++ -c src/Dimension.cpp -o Dimension.o

ScafDpData.o: src/ScafDpData.cpp
	g++ -c src/ScafDpData.cpp -o ScafDpData.o

tinymt32.o: src/TinyMT-src-1.0.1/tinymt/tinymt32.c
	g++ -c src/TinyMT-src-1.0.1/tinymt/tinymt32.c -o tinymt32.o

misc.o: src/misc.cpp
	g++ -c src/misc.cpp -o misc.o
