#
#	Makefile	written by Paul Harrison, October 1996
#

pplb:	pplb.o ppff.o
	cc -ansi -o pplb pplb.o ppff.o -lm -lnagc

pplb.o:	pplb.c
	cc -ansi -I ~/header -c pplb.c

ppff.o:	ppff.c
	cc -ansi -I ~/header -c ppff.c

ppsl:	ppsl.o ppff.o
	cc -ansi -o ppsl ppsl.o ppff.o -lm -lnagc

ppsl.o:	ppsl.c
	cc -ansi -I ~/header -c ppsl.c

ppvfq:	ppvfq.o ppff.o
	cc -ansi -o ppvfq ppvfq.o ppff.o -lm

ppvfq.o:ppvfq.c
	cc -ansi -I ~/header -c ppvfq.c

ppcd:	ppcd.c
	cc -ansi -I ~/header -o ppcd ppcd.c -lm
