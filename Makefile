APP=libgpr
VERSION=1.03
RELEASE=1
SONAME=${APP}.so.0
LIBNAME=${APP}-${VERSION}.so.0.0.${RELEASE}
ARCH_TYPE=`uname -m`

all:
	gcc -shared -Wl,-soname,${SONAME} -std=c99 -pedantic -fPIC -O3 -o ${LIBNAME} src/*.c -Isrc -lm -lz -fopenmp
debug:
	gcc -shared -Wl,-soname,${SONAME} -std=c99 -pedantic -fPIC -g -o ${LIBNAME} src/*.c -Isrc -lm -lz -fopenmp
tests:
	gcc -Wall -std=c99 -pedantic -g -o $(APP)_tests unittests/*.c src/*.c -Isrc -Iunittests -lm -lz -fopenmp

ltest:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest/*.c -lgpr -lm -lz -fopenmp

ltestc:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_cartesian/*.c -lgpr -lm -lz -fopenmp

ltestm:
	gcc -Wall -std=c99 -pedantic -g -o $(APP) libtest_morph/*.c -lgpr -lm -lz -fopenmp
source:
	tar -cvzf ../${APP}_${VERSION}.orig.tar.gz ../${APP}-${VERSION} --exclude-vcs
install:
	mkdir -p ${DESTDIR}/usr
	mkdir -p ${DESTDIR}/usr/lib
	mkdir -p ${DESTDIR}/usr/lib/${APP}
	mkdir -p ${DESTDIR}/usr/include
	mkdir -p ${DESTDIR}/usr/include/${APP}
	cp src/*.h ${DESTDIR}/usr/include/${APP}
	install -m 755 ${LIBNAME} ${DESTDIR}/usr/lib
	ln -sf ${DESTDIR}/usr/lib/${LIBNAME} ${DESTDIR}/usr/lib/${SONAME}
	ln -sf ${DESTDIR}/usr/lib/${LIBNAME} ${DESTDIR}/usr/lib/${APP}.so
	ldconfig
	mkdir -m 755 -p ${DESTDIR}/usr/share
	mkdir -m 755 -p ${DESTDIR}/usr/share/man
	mkdir -m 755 -p ${DESTDIR}/usr/share/man/man1
	install -m 644 man/${APP}.1.gz ${DESTDIR}/usr/share/man/man1
uninstall:
	rm -f /usr/share/man/man1/${APP}.1.gz
	rm -f /usr/lib/${LIBNAME}
	rm -f /usr/lib/${APP}.so
	rm -f /usr/lib/${SONAME}
	rm -rf /usr/include/${APP}
	ldconfig
instlib:
	mkdir -p ${DESTDIR}/usr
	mkdir -p ${DESTDIR}/usr/lib
	mkdir -p ${DESTDIR}/usr/lib/${APP}
	mkdir -p ${DESTDIR}/usr/include
	mkdir -p ${DESTDIR}/usr/include/${APP}
	cp src/*.h ${DESTDIR}/usr/include/${APP}
	install -m 755 ${LIBNAME} ${DESTDIR}/usr/lib
	mkdir -m 755 -p ${DESTDIR}/usr/share
	mkdir -m 755 -p ${DESTDIR}/usr/share/man
	mkdir -m 755 -p ${DESTDIR}/usr/share/man/man1
	install -m 644 man/${APP}.1.gz ${DESTDIR}/usr/share/man/man1
clean:
	rm -f ${LIBNAME} \#* \.#* gnuplot* *.png debian/*.substvars debian/*.log
	rm -fr deb.* debian/${APP} rpmpackage/${ARCH_TYPE}
	rm -f ../${APP}*.deb ../${APP}*.changes ../${APP}*.asc ../${APP}*.dsc
	rm -f rpmpackage/*.src.rpm archpackage/*.gz archpackage/*.xz
	rm -f  puppypackage/*.gz puppypackage/*.pet slackpackage/*.txz
