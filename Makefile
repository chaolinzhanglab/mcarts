# uncomment mv commands in "check.txt" section if you have write privileges to /usr/local/bin

.PHONY: all clean

all: check.txt

PatternMatch/PatternMatch:
	cd PatternMatch_src; \
	make; \
	cd ..

RegExpMatch/RegExpMatch:
	cd RegExpMatch_src; \
	make; \
	cd ..

check.txt: PatternMatch/PatternMatch RegExpMatch/RegExpMatch
	chmod 755 *.pl; \
	chmod 755 PatternMatch_src/PatternMatch; \
	chmod 755 RegExpMatch_src/RegExpMatch; \
	ln -s PatternMatch_src/PatternMatch; \
	ln -s RegExpMatch_src/RegExpMatch; \
	touch check.txt; \
	# mv PatternMatch/PatternMatch /usr/local/bin; \
	# mv RegExpMatch/RegExpMatch /usr/local/bin; \

clean:
	rm check.txt; #cd PatternMatch_src; make clean; cd ../RegExpMatch_src; make clean; cd ..

