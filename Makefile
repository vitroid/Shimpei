all: 1h.nx3a 1h.dope10.+1.+10.pickle 1h.dope10.yap 1h.dope10.mdva
%.dope10.pickle: 1h.nx3a dope.py
	python test_dope.py $< 10 $@
%.mdva: %.pickle to_mdview.py
	python to_mdview.py $< > $@
%.yap: %.pickle to_yaplot.py
	python to_yaplot.py $< > $@
%.+1.pickle: %.pickle trial_move.py
	python trial_move.py $< 1 $@
%.+10.pickle: %.pickle trial_move.py
	python trial_move.py $< 10 $@


1h.nx3a:
	genice 1h -r 4 4 4 -f e > $@

clean:
	-rm 1h.* rm *.mdva *.yap *.pickle
