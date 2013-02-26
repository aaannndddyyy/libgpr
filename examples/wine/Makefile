
all:
	gcc -Wall -ansi -pedantic -o winequality winequality.c -lgpr -lm -lz -fopenmp

debug:
	gcc -Wall -ansi -pedantic -g -o winequality winequality.c -lgpr -lm -lz -fopenmp

clean:
	rm -f *.o winequality agent agent.c server.rb
