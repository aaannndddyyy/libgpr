
all:
	gcc -Wall -ansi -pedantic -o concreteslump concreteslump.c -lgpr -lm -lz -fopenmp

debug:
	gcc -Wall -ansi -pedantic -g -o concreteslump concreteslump.c -lgpr -lm -lz -fopenmp

clean:
	rm -f *.o concreteslump agent agent.c server.rb *.png *.dot
