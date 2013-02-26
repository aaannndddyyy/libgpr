
all:
	gcc -Wall -ansi -pedantic -o parkinsons parkinsons.c -lgpr -lm -lz -fopenmp

debug:
	gcc -Wall -ansi -pedantic -g -o parkinsons parkinsons.c -lgpr -lm -lz -fopenmp

clean:
	rm -f *.o parkinsons agent agent.c server.rb
