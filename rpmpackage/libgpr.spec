Name: libgpr
Version: 1.03
Release: 1%{?dist}
Summary: Library for genetic programming
License: BSD
URL: https://github.com/fuzzgun/libgpr
Packager: Bob Mottram (4096 bits) <bob@robotics.uk.to>
Source0: http://yourdomainname.com/src/%{name}_%{version}.orig.tar.gz
Group: Development/ArtificialIntelligence

Requires: gnuplot, zlib


%description
Making the inclusion of Genetic Programming easy within any C/C++
application. Genetic programming (GP) is a powerful technique, inspired by
the process of natural selection, which can be utilized to automatically
discover programs which produce a desired input to output transformation.
Both classical tree-based and Cartesian forms of Genetic Programming are
supported, including self-modifying variants.

%prep
%setup -q

%build
%configure
make %{?_smp_mflags}

%install
rm -rf %{buildroot}
mkdir -p %{buildroot}
mkdir -p %{buildroot}/etc
mkdir -p %{buildroot}/etc/%{name}
mkdir -p %{buildroot}/usr
mkdir -p %{buildroot}/usr/bin
mkdir -p %{buildroot}/usr/lib
mkdir -p %{buildroot}/usr/lib/%{name}
mkdir -p %{buildroot}/usr/share
mkdir -p %{buildroot}/usr/share/man
mkdir -p %{buildroot}/usr/share/man/man1
# Make install but to the RPM BUILDROOT directory
make instlib -B DESTDIR=%{buildroot}

%files
%doc README.md LICENSE
%defattr(-,root,root,-)
%{_bindir}/*
%{_mandir}/man1/*

%post
umask 007
ldconfig > /dev/null 2>&1

%postun
umask 007
ldconfig > /dev/null 2>&1

%changelog
* Thu Mar 21 2013 Bob Mottram (4096 bits) <bob@robotics.uk.to> - 1.03-1
- Programs can now use complex numbers

* Thu Mar 7 2013 Bob Mottram (4096 bits) <bob@robotics.uk.to> - 1.02-1
- Added chromosomes to Cartesian genetic programming
- Example applications
- Dropouts to reduce overfitting
- Custom functions can be defined
- Added environment structure
- Added artificial life example
- Added leaf classification example
- Added genetic art example

* Mon Nov 19 2012 Bob Mottram (4096 bits) <bob@robotics.uk.to> - 1.01-1
- Bug fix for fitness history graph
- Addition of Automatically Defined Functions

* Sun Nov 18 2012 Bob Mottram (4096 bits) <bob@robotics.uk.to> - 1.00-1
- Initial release
