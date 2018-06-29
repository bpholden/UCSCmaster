SHELL = /bin/sh
INSTALL = install -c
MKDIR = install -d 

destdir = ../bin/

srcdir = .

SRC =  	make_googledexdat.py NightSim.py sim_summary.py update_googledex.py pri_restore.py sim_night.py sim_nights.py \
	calc_cadence.py calc_etime_precision.py calc_precision_mag.py conv_bj.py download_googledex.py exposure_time_test.py \
	Generate_Errors.py gen_starlist_entry.py gen_template_entry.py make_owner_list.py make_scriptobsline.py \
	make_template_list.py run_master.py master.config



install:
	$(MKDIR) $(destdir)
	$(INSTALL) $(SRC) $(destdir)
