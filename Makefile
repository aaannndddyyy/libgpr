APP=libgpr
VERSION=1.03
RELEASE=1
SONAME=${APP}.so.0
LIBNAME=${APP}-${VERSION}.so.0.0.${RELEASE}
ARCH_TYPE=`uname -m`
PREFIX?=/usr/local

all:
	gcc -shared -Wl,-soname,${SONAME} -std=c99 -pedantic -fPIC -O3 -o ${LIBNAME} src/*.c -Isrc -lm -lz -fopenmp
debug:
	gcc -shared -Wl,-soname,${SONAME} -std=c99 -pedantic -fPIC -g -o ${LIBNAME} src/*.c -Isrc -lm -lz -fopenmp
source:
	tar -cvf ../${APP}_${VERSION}.orig.tar ../${APP}-${VERSION} --exclude-vcs
	gzip -f9n ../${APP}_${VERSION}.orig.tar
install:
	mkdir -p ${DESTDIR}/usr
	mkdir -p ${DESTDIR}${PREFIX}
	mkdir -p ${DESTDIR}${PREFIX}/lib
	mkdir -p ${DESTDIR}${PREFIX}/lib/${APP}
	mkdir -p ${DESTDIR}${PREFIX}/include
	mkdir -p ${DESTDIR}${PREFIX}/include/${APP}
	cp src/*.h ${DESTDIR}${PREFIX}/include/${APP}
	install -m 755 ${LIBNAME} ${DESTDIR}${PREFIX}/lib
	ln -sf ${DESTDIR}${PREFIX}/lib/${LIBNAME} ${DESTDIR}${PREFIX}/lib/${SONAME}
	ln -sf ${DESTDIR}${PREFIX}/lib/${LIBNAME} ${DESTDIR}${PREFIX}/lib/${APP}.so
	ldconfig
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share/man
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share/man/man1
	install -m 644 man/${APP}.1.gz ${DESTDIR}${PREFIX}/share/man/man1
uninstall:
	rm -f ${PREFIX}/share/man/man1/${APP}.1.gz
	rm -f ${PREFIX}/lib/${LIBNAME}
	rm -f ${PREFIX}/lib/${APP}.so
	rm -f ${PREFIX}/lib/${SONAME}
	rm -rf ${PREFIX}/include/${APP}
	ldconfig
instlib:
	mkdir -p ${DESTDIR}/usr
	mkdir -p ${DESTDIR}${PREFIX}
	mkdir -p ${DESTDIR}${PREFIX}/lib
	mkdir -p ${DESTDIR}${PREFIX}/lib/${APP}
	mkdir -p ${DESTDIR}${PREFIX}/include
	mkdir -p ${DESTDIR}${PREFIX}/include/${APP}
	cp src/*.h ${DESTDIR}${PREFIX}/include/${APP}
	install -m 755 ${LIBNAME} ${DESTDIR}${PREFIX}/lib
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share/man
	mkdir -m 755 -p ${DESTDIR}${PREFIX}/share/man/man1
	install -m 644 man/${APP}.1.gz ${DESTDIR}${PREFIX}/share/man/man1
clean:
	rm -f ${LIBNAME} \#* \.#* gnuplot* *.png debian/*.substvars debian/*.log
	rm -fr deb.* debian/${APP} rpmpackage/${ARCH_TYPE}
	rm -f ../${APP}*.deb ../${APP}*.changes ../${APP}*.asc ../${APP}*.dsc
	rm -f rpmpackage/*.src.rpm archpackage/*.gz archpackage/*.xz
	rm -f puppypackage/*.gz puppypackage/*.pet slackpackage/*.txz

sourcedeb:
	tar -cvf ../${APP}_${VERSION}.orig.tar ../${APP}-${VERSION} --exclude-vcs --exclude 'debian'
	gzip -f9n ../${APP}_${VERSION}.orig.tar
