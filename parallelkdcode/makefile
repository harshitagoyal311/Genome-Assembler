exec: main.o parallelkdtree.o kdtree.o utils.o data.o breakfile
	mpicc -o main main.o parallelkdtree.o kdtree.o utils.o data.o -lm

breakfile: 
	gcc breakfile.c data.c -o breakfile -lm

main.o: main.c
	mpicc -c main.c -lm

parallelkdtree.o: parallelkdtree.c
	mpicc -c parallelkdtree.c -lm

kdtree.o: kdtree.c
	mpicc -c kdtree.c -lm

utils.o: utils.c
	mpicc -c utils.c -lm

data.o: data.c
	mpicc -c data.c -lm

clean:
	rm -f *.o
	rm -f out*
	rm -f test*
	rm -f breakfile