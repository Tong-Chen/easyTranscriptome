# easyTranscriptome

Utils.

CPC2: https://github.com/gao-lab/CPC2_standalone/releases/tag/v1.0.1
	* With line 311 in CPC2.py changing 
	from
		script_dir,filename = os.path.split(os.path.abspath(sys.argv[0]))
	to 
		script_dir,filename = os.path.split(os.path.realpath(sys.argv[0]))

