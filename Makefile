all: configure configure.win tools/config.rpath

tools/config.rpath:
	mkdir -p tools
	touch tools/config.rpath

configure: configure.ac tools/config.rpath 
	autoreconf -fi -I tools
	rm -f aclocal.m4 

configure.win: configure
	cp configure configure.win

clean:
	rm -f src/Makevars src/Makevars.ucrt 

.PHONY: all