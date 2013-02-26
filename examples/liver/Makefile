
all:
	gcc -Wall -ansi -pedantic -o liver liver.c -lgpr -lm -lz -fopenmp

debug:
	gcc -Wall -ansi -pedantic -g -o liver liver.c -lgpr -lm -lz -fopenmp

clean:
	rm -f *.o liver agent agent.c server.rb
