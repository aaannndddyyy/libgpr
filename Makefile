APP=libgpr
VERSION=1.02
RELEASE=1
SONAME=$(APP).so.0
LIBNAME=$(APP)-$(VERSION).so.0.0.$(RELEASE)
USRBASE=/usr

all:
	gcc -c -std=c99 -pedantic -fPIC -o som.o src/som.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP).o src/gpr.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP)c.o src/gprc.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP)cm.o src/gprcm.c -Isrc -lm -fopenmp
	gcc -shared -Wl,-soname,$(SONAME) -o $(LIBNAME) som.o $(APP).o $(APP)c.o $(APP)cm.o
#	objdump -p ${LIBNAME} | sed -n -e's/^[[:space:]]*SONAME[[:space:]]*//p' | sed -e's/\([0-9]\)\.so\./\1-/; s/\.so\.//'

debug:
	gcc -c -std=c99 -pedantic -fPIC -g -o som.o src/som.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP).o src/gpr.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP)c.o src/gprc.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP)cm.o src/gprcm.c -Isrc -lm -fopenmp
	gcc -shared -Wl,-soname,$(SONAME) -o $(LIBNAME) som.o $(APP).o $(APP)c.o $(APP)cm.o

tests:
	gcc -Wall -std=c99 -pedantic -g -o $(APP)_tests unittests/*.c src/*.c -Isrc -Iunittests -lm -fopenmp

ltest:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest/*.c -lgpr -lm -fopenmp

ltestc:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_cartesian/*.c -lgpr -lm -fopenmp

ltestm:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_morph/*.c -lgpr -lm -fopenmp

source:
	tar -cvzf ../$(APP)_$(VERSION).orig.tar.gz ../$(APP)-$(VERSION) --exclude=.bzr

install:
	mkdir -m 755 -p $(USRBASE)/lib
	cp $(LIBNAME) $(USRBASE)/lib
	cp man/$(APP).1.gz $(USRBASE)/share/man/man1
	mkdir -m 755 -p $(USRBASE)/include/$(APP)
	cp src/*.h $(USRBASE)/include/$(APP)
	chmod 755 $(USRBASE)/lib/$(LIBNAME)
	chmod 644 $(USRBASE)/include/$(APP)/*.h
	chmod 644 $(USRBASE)/share/man/man1/$(APP).1.gz
	ln -sf $(USRBASE)/lib/$(LIBNAME) $(USRBASE)/lib/$(SONAME)
	ln -sf $(USRBASE)/lib/$(LIBNAME) $(USRBASE)/lib/$(APP).so
	ldconfig

clean:
	rm -f examples/random/*.png examples/random/*.dot examples/random/random
	rm -f examples/pursuer_evader/*.png examples/pursuer_evader/*.dot examples/pursuer_evader/pursuer
	rm -f examples/artificial_ant/*.png examples/artificial_ant/*.dot examples/artificial_ant/ant
	rm -f examples/cart_centering/*.png examples/cart_centering/*.dot examples/cart_centering/cart
	rm -f examples/economic_modeling/*.png examples/economic_modeling/*.dot examples/economic_modeling/econ
	rm -f examples/symbolic_regression/*.png examples/symbolic_regression/*.dot examples/symbolic_regression/symreg
	rm -f $(LIBNAME) $(APP) $(APP).so.* $(APP).o $(APP)_tests $(APP)-* *.dot
	rm -f *.dat *.png *.txt *.rb agent.c agent
	rm -f \#* \.#* debian/*.substvars debian/*.log *.so.0.0.1 *.o
	rm -rf deb.* debian/$(APP)0 debian/$(APP)0-dev
	rm -f ../$(APP)*.deb ../$(APP)*.changes ../$(APP)*.asc ../$(APP)*.dsc ../$(APP)_$(VERSION)*.gz

