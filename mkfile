MKSHELL=$PLAN9/bin/rc

<|echo 'gps ''='' *.gp'
pdfs = ${gps:%.gp=%.pdf}

all:V: $pdfs
%.pdf:	%.gp
	gpctx -s 11in,8.5in $stem.gp
