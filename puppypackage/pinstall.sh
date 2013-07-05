#!/bin/bash

APP=libgpr
VERSION=1.03
RELEASE=1
SONAME=$(APP).so.0
LIBNAME=$(APP)-$(VERSION).so.0.0.$(RELEASE)

ln -sf /usr/lib/$(LIBNAME) /usr/lib/$(SONAME)
ln -sf /usr/lib/$(LIBNAME) /usr/lib/$(APP).so
ldconfig
