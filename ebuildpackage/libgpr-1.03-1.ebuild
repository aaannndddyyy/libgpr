# $Header: $

EAPI=4

DESCRIPTION="Making the inclusion of Genetic Programming easy within any C/C++ application. Genetic programming (GP) is a powerful technique, inspired by the process of natural selection, which can be utilized to automatically discover programs which produce a desired input to output transformation. Both classical tree-based and Cartesian forms of Genetic Programming are supported, including self-modifying variants."
HOMEPAGE="https://github.com/fuzzgun/libgpr"
SRC_URI="${PN}/${P}.tar.gz"
LICENSE="BSD"
SLOT="0"
KEYWORDS="x86"
RDEPEND="dev-libs/popt"
DEPEND="${RDEPEND}"

src_configure() {
    econf --with-popt
}

src_install() {
    emake DESTDIR="${D}" instlib
    # Install README and (Debian) changelog
    dodoc README.md debian/changelog
}
