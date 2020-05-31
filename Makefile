CXX = g++

MAIN = meta_platanus.pl
SUBDIRS = src/minimap2 src/meta_platanus src/meta_platanus_phase src/megahit src/samtools src/bwa
SUBBIN = sub_bin

.PHONY:all clean $(SUBDIRS) $(SUBBIN)


all: $(MAIN) $(SUBDIRS) $(SUBBIN)

$(MAIN):

src/samtools:
	cd src/samtools && ./configure

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBBIN):
	mkdir -p $(SUBBIN)
	cp src/meta_platanus/meta_platanus $(SUBBIN)
	cp src/meta_platanus_phase/meta_platanus_phase $(SUBBIN)
	cp src/minimap2/minimap2 $(SUBBIN)
	cp src/megahit/megahit* $(SUBBIN)
	cp src/samtools/samtools $(SUBBIN)
	cp src/bwa/bwa $(SUBBIN)
	cp src/scripts/*.pl src/scripts/*.pm $(SUBBIN)
	cp pre_built_bin/metabat/metabat2 $(SUBBIN)
	cp pre_built_bin/metabat/jgi_summarize_bam_contig_depths $(SUBBIN)

clean:
	rm -rf $(SUBBIN)
	cd src/meta_platanus && make clean
	cd src/meta_platanus_phase && make clean
	cd src/minimap2 && make clean
	cd src/megahit && make clean
	cd src/samtools && make clean
	cd src/bwa&& make clean
#	$(MAKE) clean -C $(SUBDIRS)
