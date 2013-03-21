APP=libgpr
VERSION=1.02
RELEASE=1
SONAME=$(APP).so.0
LIBNAME=$(APP)-$(VERSION).so.0.0.$(RELEASE)
USRBASE=/usr

all:
	gcc -c -std=c99 -pedantic -fPIC -o pnglite.o src/pnglite.c -Isrc
	gcc -c -std=c99 -pedantic -fPIC -o som.o src/som.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP).o src/gpr.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP)c.o src/gprc.c -Isrc -lm -lz -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -o $(APP)cm.o src/gprcm.c -Isrc -lm -lz -fopenmp
	gcc -shared -Wl,-soname,$(SONAME) -o $(LIBNAME) pnglite.o som.o $(APP).o $(APP)c.o $(APP)cm.o
#	objdump -p ${LIBNAME} | sed -n -e's/^[[:space:]]*SONAME[[:space:]]*//p' | sed -e's/\([0-9]\)\.so\./\1-/; s/\.so\.//'

debug:
	gcc -c -std=c99 -pedantic -fPIC -g -o pnglite.o src/pnglite.c -Isrc
	gcc -c -std=c99 -pedantic -fPIC -g -o som.o src/som.c -Isrc -lm -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP).o src/gpr.c -Isrc -lm -lz -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP)c.o src/gprc.c -Isrc -lm -lz -fopenmp
	gcc -c -std=c99 -pedantic -fPIC -g -o $(APP)cm.o src/gprcm.c -Isrc -lm -lz -fopenmp
	gcc -shared -Wl,-soname,$(SONAME) -o $(LIBNAME) pnglite.o som.o $(APP).o $(APP)c.o $(APP)cm.o

tests:
	gcc -Wall -std=c99 -pedantic -g -o $(APP)_tests unittests/*.c src/*.c -Isrc -Iunittests -lm -lz -fopenmp

ltest:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest/*.c -lgpr -lm -fopenmp

ltestc:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_cartesian/*.c -lgpr -lm -lz -fopenmp

ltestm:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_morph/*.c -lgpr -lm -lz -fopenmp

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
	rm -f examples/arrhythmia/*.png examples/arrhythmia/*.dot examples/arrhythmia/arrhythmia examples/arrhythmia/*.rb examples/arrhythmia/agent.c examples/arrhythmia/agent
	rm -f examples/art/*.png examples/art/*.dot examples/art/art examples/art/*.rb examples/art/agent.c examples/art/agent
	rm -f examples/artificial_ant/*.png examples/artificial_ant/*.dot examples/artificial_ant/ant examples/artificial_ant/*.rb examples/artificial_ant/agent.c examples/artificial_ant/agent
	rm -f examples/cancer_classification/*.png examples/cancer_classification/*.dot examples/cancer_classification/cancer.c examples/cancer_classification/*.rb examples/cancer_classification/agent.c examples/cancer_classification/agent
	rm -f examples/cart_centering/*.png examples/cart_centering/*.dot examples/cart_centering/cart examples/cart_centering/*.rb examples/cart_centering/agent.c examples/cart_centering/agent
	rm -f examples/concreteslump/*.png examples/concreteslump/*.dot examples/concreteslump/concreteslump examples/concreteslump/*.rb examples/concreteslump/agent.c examples/concreteslump/agent
	rm -f examples/critters/*.png examples/critters/*.dot examples/critters/critters examples/critters/agent.c examples/critters/agent
	rm -f examples/leaves/*.png examples/leaves/*.dot examples/leaves/leaves examples/leaves/species* examples/leaves/*.rb examples/leaves/agent.c examples/leaves/agent
	rm -f examples/liver/*.png examples/liver/*.dot examples/liver/liver examples/liver/*.rb examples/liver/agent.c examples/liver/agent
	rm -f examples/parkinsons/*.png examples/parkinsons/*.dot examples/parkinsons/parkinsons examples/parkinsons/*.rb examples/parkinsons/agent.c examples/parkinsons/agent
	rm -f examples/pursuer_evader/*.png examples/pursuer_evader/*.dot examples/pursuer_evader/pursuer examples/pursuer_evader/*.rb examples/pursuer_evader/agent.c examples/pursuer_evader/agent
	rm -f examples/random/*.png examples/random/*.dot examples/random/random examples/random/*.rb examples/random/agent.c examples/random/agent
	rm -f examples/symbolic_regression/*.png examples/symbolic_regression/*.dot examples/symbolic_regression/symreg examples/symbolic_regression/*.rb examples/symbolic_regression/agent.c examples/symbolic_regression/agent
	rm -f examples/economic_modeling/*.png examples/economic_modeling/*.dot examples/economic_modeling/econ examples/economic_modeling/*.rb examples/economic_modeling/agent.c examples/economic_modeling/agent
	rm -f examples/violentcrime/*.png examples/violentcrime/*.dot examples/violentcrime/violentcrime examples/violentcrime/*.rb examples/violentcrime/agent.c examples/violentcrime/agent
	rm -f examples/wine/*.png examples/wine/*.dot examples/wine/wine examples/wine/*.rb examples/wine/agent.c examples/wine/agent
	rm -f $(LIBNAME) $(APP) $(APP).so.* $(APP).o $(APP)_tests $(APP)-* *.dot
	rm -f *.dat *.png *.txt *.rb agent.c agent
	rm -f \#* \.#* debian/*.substvars debian/*.log *.so.0.0.1 *.o temp_agent
	rm -rf deb.* debian/$(APP)0 debian/$(APP)0-dev
	rm -f ../$(APP)*.deb ../$(APP)*.changes ../$(APP)*.asc ../$(APP)*.dsc ../$(APP)_$(VERSION)*.gz

