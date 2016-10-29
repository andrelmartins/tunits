all:
	rm -Rf pkg.roxygen
	(R CMD roxygen pkg)||(./Roxygen pkg)
	R CMD INSTALL pkg.roxygen

clean:
	rm -Rf pkg.roxygen
