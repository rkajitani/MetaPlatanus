CXX = g++

MAIN = metaplatanus
SUBDIRS = src/metaplatanus_core
SUBBIN = sub_bin

.PHONY:all clean $(SUBDIRS) $(SUBBIN)


all: $(MAIN) $(SUBDIRS) $(SUBBIN)

$(MAIN):

$(SUBDIRS):
	$(MAKE) -C $@

$(SUBBIN):
	mkdir -p $(SUBBIN)
	cp src/metaplatanus_core/metaplatanus_core $(SUBBIN)
	cp src/scripts/*.pl $(SUBBIN)

clean:
	rm -rf $(SUBBIN)
	cd src/metaplatanus_core && make clean
#	$(MAKE) clean -C $(SUBDIRS)
