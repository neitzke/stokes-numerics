# -*-mode: Makefile-*-

PYINTERP=python3
SCRIPTDIR=
OUTDIR=output

PFH=$(PYINTERP) $(SCRIPTDIR)paperfig-hitchin.py
PFO=$(PYINTERP) $(SCRIPTDIR)paperfig-oper.py
PFHEAT=$(PYINTERP) $(SCRIPTDIR)paperfig-heatmaps.py
PFHK=$(PYINTERP) $(SCRIPTDIR)paperfig-hkmetric.py
PFXH=$(PYINTERP) $(SCRIPTDIR)paperfig-xarshitchin.py
PFXO=$(PYINTERP) $(SCRIPTDIR)paperfig-xaroper.py

PTO=$(PYINTERP) $(SCRIPTDIR)papertable-oper.py
PTH=$(PYINTERP) $(SCRIPTDIR)papertable-hitchin.py

HITCHINPDFS=\
	$(OUTDIR)/paperfigA1A2th0.0small.pdf \
	$(OUTDIR)/paperfigA1A2th0.0large.pdf \
	$(OUTDIR)/paperfigA1A2th0.1small.pdf \
	$(OUTDIR)/paperfigA1A2th0.1large.pdf \
	$(OUTDIR)/paperfigA1A2th0.0-deformedsmall.pdf \
	$(OUTDIR)/paperfigA1A2th0.0-deformedlarge.pdf \
	$(OUTDIR)/paperfigA1A3th0.1small.pdf \
	$(OUTDIR)/paperfigA1A3th0.1large.pdf \
	$(OUTDIR)/paperfigA2A1th0.0small.pdf \
	$(OUTDIR)/paperfigA2A1th0.0large.pdf \
	$(OUTDIR)/paperfigA2A2th0.1small.pdf \
	$(OUTDIR)/paperfigA2A2th0.1large.pdf

OPERPDFS=\
	$(OUTDIR)/paperfigA1A2operth0.0large.pdf \
	$(OUTDIR)/paperfigA1A2operth0.0small.pdf \
	$(OUTDIR)/paperfigA1A2operth0.1large.pdf \
	$(OUTDIR)/paperfigA1A2operth0.1small.pdf \
	$(OUTDIR)/paperfigA1A2operth0.0-deformedsmall.pdf \
	$(OUTDIR)/paperfigA1A2operth0.0-deformedlarge.pdf \
	$(OUTDIR)/paperfigA1A3operth0.1small.pdf \
	$(OUTDIR)/paperfigA1A3operth0.1large.pdf \
	$(OUTDIR)/paperfigA2A1operth0.0small.pdf \
	$(OUTDIR)/paperfigA2A1operth0.0large.pdf \
	$(OUTDIR)/paperfigA2A1deformedoperth0.0small.pdf \
	$(OUTDIR)/paperfigA2A1deformedoperth0.0large.pdf \
	$(OUTDIR)/paperfigA2A2operth0.1small.pdf \
	$(OUTDIR)/paperfigA2A2operth0.1large.pdf

XARPDFS=$(OUTDIR)/paperfigxaroper.pdf $(OUTDIR)/paperfigxarshitchin.pdf

HKPDFS=$(OUTDIR)/paperfigheatmaps.pdf $(OUTDIR)/paperfighkmetric.pdf

PDFS=$(HITCHINPDFS) $(OPERPDFS) $(XARPDFS) $(HKPDFS) 

TABLES=$(OUTDIR)/table-a1a2deformed-opers.inc $(OUTDIR)/table-a1a2-hitchin.inc

# Options applied to all hitchin and oper plots
COMMONOPTS=--reldiff --diffymin=5e-16 --diffymax=5

# Extra options for all hitchin plots
HITCHINOPTS=--pde-nmesh 2047,4095,8191

# Extra options for all oper plots
OPEROPTS=

all: hitchin oper hk xar tables

$(OUTDIR):
	@[ -d $@ ] || ( echo "mkdir -p $@"; mkdir -p $@ )

hitchin: $(HITCHINPDFS)

oper: $(OPERPDFS)

hk: $(HKPDFS)

xar: $(XARPDFS)

tables: $(TABLES)

clean:
	rm -f $(PDFS) $(TABLES)

clean-hitchin:
	rm -f $(HITCHINPDFS)

clean-oper:
	rm -f $(OPERPDFS)

clean-hk:
	rm -f $(HKPDFS)

clean-tables:
	rm -f $(TABLES)

clean-xar:
	rm -f $(XARPDFS)

.PHONY: all hitchin oper hk tables clean clean-hitchin clean-oper clean-hk clean-xar clean-tables $(OUTDIR)

# ----------------------------------------------------------------------
#                       HITCHIN SECTION PLOTS
# ----------------------------------------------------------------------

$(OUTDIR)/paperfigA1A2th0.0small.pdf: | $(OUTDIR)
	$(PFH) A1A2 0.0 --xlog --ylinear --rmax 0.1 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2th0.0large.pdf: | $(OUTDIR)
	$(PFH) A1A2 0.0 --xlinear --ylog --rmin 0.1 --rmax 10 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2th0.1small.pdf: | $(OUTDIR)
	$(PFH) A1A2 0.1 --xlog --ylinear --rmax 0.1 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2th0.1large.pdf: | $(OUTDIR)
	$(PFH) A1A2 0.1 --xlinear --ylog --rmin 0.1 --rmax 10 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2th0.0-deformedsmall.pdf: | $(OUTDIR)
	$(PFH) 'A1A2_Lambda=0.8j_cr=1_ci=0' 0.0 --tex-theoryname '$$(A_1,A_2)$$ at $$\Lambda=0.8i$$, $$c=1$$, $$\vartheta=0$$' --xlog --ylinear --rmax 0.1 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2th0.0-deformedlarge.pdf: | $(OUTDIR)
	$(PFH) 'A1A2_Lambda=0.8j_cr=1_ci=0' 0.0 --tex-theoryname '$$(A_1,A_2)$$ at $$\Lambda=0.8i$$, $$c=1$$, $$\vartheta=0$$' --xlinear --ylog --rmin 0.1 --rmax 10 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A3th0.1small.pdf: | $(OUTDIR)
	$(PFH) A1A3 0.1 --exclude X3 --xlog --ylinear --rmax 0.1 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A3th0.1large.pdf: | $(OUTDIR)
	$(PFH) A1A3 0.1 --exclude X3 --xlinear --ylog --rmin 0.1 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1th0.0small.pdf: | $(OUTDIR)
	$(PFH) A2A1 0.0 --xlog --ylinear --rmax 0.1 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1th0.0large.pdf: | $(OUTDIR)
	$(PFH) A2A1 0.0 --xlinear --ylog --rmin 0.1 --rmax 8 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A2th0.1small.pdf: | $(OUTDIR)
	$(PFH) A2A2 0.1 --legend-ncol 2 --xlog --ylinear --rmax 0.1 --show-semiflat X3 X4 --title-suffix ' for small $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A2th0.1large.pdf: | $(OUTDIR)
	$(PFH) A2A2 0.1 --extra-left-yrange 0.18 --legend-ncol 2 --xlinear --ylog --rmin 0.1 --rmax 3 --show-semiflat all --title-suffix ' for large $$R$$' --output $@ $(HITCHINOPTS) $(COMMONOPTS)


# ----------------------------------------------------------------------
#                       OPER PLOTS
# ----------------------------------------------------------------------

$(OUTDIR)/paperfigA1A2operth0.0small.pdf: | $(OUTDIR)
	$(PFO) A1A2 0.0 --hmin 10 --xlog --ylinear --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2operth0.0large.pdf: | $(OUTDIR)
	$(PFO) A1A2 0.0 --hmax 10 --xlinear --ylog --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2operth0.1small.pdf: | $(OUTDIR)
	$(PFO) A1A2 0.1 --hmin 10 --xlog --ylinear --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2operth0.1large.pdf: | $(OUTDIR)
	$(PFO) A1A2 0.1 --hmax 10 --xlinear --ylog --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2operth0.0-deformedsmall.pdf: | $(OUTDIR)
	$(PFO) 'A1A2_Lambda=0.8j_cr=1_ci=0' 0.0 --tex-theoryname '$$(A_1,A_2)$$ opers at $$\Lambda=0.8i$$, $$c=1$$, $$\vartheta=0$$' --xlog --ylinear --hmin 10 --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A2operth0.0-deformedlarge.pdf: | $(OUTDIR)
	$(PFO) 'A1A2_Lambda=0.8j_cr=1_ci=0' 0.0 --tex-theoryname '$$(A_1,A_2)$$ opers at $$\Lambda=0.8i$$, $$c=1$$, $$\vartheta=0$$' --xlinear --ylog --hmax 10 --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A3operth0.1small.pdf: | $(OUTDIR)
	$(PFO) A1A3 0.1 --exclude X3 --xlog --ylinear --hmin 1 --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA1A3operth0.1large.pdf: | $(OUTDIR)
	$(PFO) A1A3 0.1 --exclude X3 --xlinear --ylog --hmax 1 --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1operth0.0small.pdf: | $(OUTDIR)
	$(PFO) A2A1 0.0 --xlog --ylinear --hmin 1 --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1operth0.0large.pdf: | $(OUTDIR)
	$(PFO) A2A1 0.0 --xlinear --ylog --hmax 1 --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1deformedoperth0.0small.pdf: | $(OUTDIR)
	$(PFO) A2A1_c=0.5 0.0 --xlog --ylinear --hmin 1 --title-suffix ' for small $$\hbar^{-1}$$' --tex-theoryname '$$(A_2,A_1)$$ opers at $$c=0.5$$, $$\vartheta=0$$' --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A1deformedoperth0.0large.pdf: | $(OUTDIR)
	$(PFO) A2A1_c=0.5 0.0 --xlinear --ylog --hmax 1 --hmin 0.1 --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --tex-theoryname '$$(A_2,A_1)$$ opers at $$c=0.5$$, $$\vartheta=0$$' --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A2operth0.1small.pdf: | $(OUTDIR)
	$(PFO) A2A2 0.1 --legend-ncol 2 --xlog --ylinear --hmin 10 --show-wkb X3 X4 --title-suffix ' for small $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

$(OUTDIR)/paperfigA2A2operth0.1large.pdf: | $(OUTDIR)
	$(PFO) A2A2 0.1 --extra-left-yrange 0.18 --legend-ncol 2 --xlinear --ylog --hmax 10 --hmin 0.3 --show-wkb all --title-suffix ' for large $$\hbar^{-1}$$' --output $@ $(OPEROPTS) $(COMMONOPTS)

# ----------------------------------------------------------------------
#                       HKMETRIC PLOTS
# ----------------------------------------------------------------------

$(OUTDIR)/paperfigheatmaps.pdf: | $(OUTDIR)
	$(PFHEAT) --multiplerows --pde-nmesh 1000 --trimfactor 0.17 --output $@

$(OUTDIR)/paperfighkmetric.pdf: | $(OUTDIR)
	$(PFHK) --Lambda 0.0 --output $@

# ----------------------------------------------------------------------
#                        XAR PLOTS
# ----------------------------------------------------------------------
$(OUTDIR)/paperfigxarshitchin.pdf: | $(OUTDIR)
	$(PFXH) --output $@

$(OUTDIR)/paperfigxaroper.pdf: | $(OUTDIR)
	$(PFXO) --output $@

# ----------------------------------------------------------------------
#                        TABLES
# ----------------------------------------------------------------------

$(OUTDIR)/table-a1a2-hitchin.inc: | $(OUTDIR)
	$(PTH) > $@

$(OUTDIR)/table-a1a2deformed-opers.inc: | $(OUTDIR)
	$(PTO) > $@
