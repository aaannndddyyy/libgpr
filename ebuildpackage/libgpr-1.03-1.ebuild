# $Header: $

EAPI=5

inherit git-2

DESCRIPTION="Making the inclusion of Genetic Programming easy within any C/C++ application. Genetic programming (GP) is a powerful technique, inspired by the process of natural selection, which can be utilized to automatically discover programs which produce a desired input to output transformation. Both classical tree-based and Cartesian forms of Genetic Programming are supported, including self-modifying variants."
HOMEPAGE="https://github.com/bashrc/libgpr"
EGIT_REPO_URI="https://github.com/bashrc/libgpr.git"
LICENSE="BSD"
SLOT="0"
KEYWORDS="x86"
DEPEND="dev-libs/popt"
RDEPEND="${DEPEND}"

src_configure() {
    econf --with-popt
}

src_compile() { :; }

src_install() {
    emake DESTDIR="${D}" PREFIX="/usr" instlib
    # Install README and (Debian) changelog
    dodoc README.md debian/changelog
}
