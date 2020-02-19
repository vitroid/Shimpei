%.dope10.pickle: 1h.nx3a dope.py
	python dope.py $< 10 $@
%.mdva: %.pickle
	python to_mdview.py $< > $@
%.+1.pickle: %.pickle
	python trial_move.py $< 1 $@


1h.nx3a:
	genice 1h -r 4 4 4 -f e > $@
