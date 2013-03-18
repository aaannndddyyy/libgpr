#
# rpm spec for libgpr
#

%define        __spec_install_post %{nil}
%define          debug_package %{nil}
%define        __os_install_post %{_dbpath}/brp-compress

Summary: Genetic programming
Name: libgpr
Version: 1.02
Release: 1
License: BSD
Group: libs
SOURCE0 : %{name}-%{version}.tar.gz
URL: https://launchpad.net/libgpr
Packager: Bob Mottram <bob@sluggish.dyndns.org>
Requires: gnuplot
BuildRoot: %{_tmppath}/%{name}-%{version}-%{release}-root

%description
The aim of libgpr is to make using Genetic Programming easy to include within any C/C++ application. Genetic programming (GP) is a powerful technique, inspired by the process of natural selection, which can be utilized to automatically discover programs which produce a desired input to output transformation. The problem can be stated simply as:

   "Given some behavior find a program which produces this behavior, or a close approximation of it"

Both classical tree-based and Cartesian forms of Genetic Programming are supported, including self-modifying CGP.

%prep
%setup -q

%build
 # Empty section.

%install
rm -rf %{buildroot}
mkdir -p  %{buildroot}

# in builddir
cp -a * %{buildroot}

%clean
rm -rf %{buildroot}

%files
%defattr(-,root,root,-)
%config(noreplace) %{_sysconfdir}/%{name}/%{name}.conf
%attr(644,root,root) /usr/lib/*
%attr(644,root,root) /usr/include/%{name}/*
%attr(644,root,root) /usr/share/man/man1/%{name}.1.gz

%post
umask 007
ldconfig > /dev/null 2>&1
ln -sf /usr/lib/%{name}-%{version}.so.0.0.1 /usr/lib/%{name}.so

%postun
umask 007
ldconfig > /dev/null 2>&1
rm /usr/lib/%{name}.so

%changelog
* Thu Nov 8 2012  Bob Mottram <bob@sluggish.dyndns.org>
- Spec file created

