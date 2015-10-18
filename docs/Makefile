all:
	latex fastcluster.Rtex
	latex fastcluster.Rtex
	latex fastcluster.Rtex
	dvipdfmx fastcluster.dvi
	mkdir keep
	mv fastcluster.pdf keep
	mv fastcluster.Rtex keep
	rm fastcluster.*
	mv keep/* .
	rmdir keep
