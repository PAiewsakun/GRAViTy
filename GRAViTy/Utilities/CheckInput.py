import optparse, os

def check_FILEPATH (option, opt_str, value, parser):
	if value != None and (not os.path.isfile(value)):
		raise optparse.OptionValueError("%s doesn't exist."%value)
	setattr(parser.values, option.dest, value)

def check_FILEPATHS (option, opt_str, values, parser):
	if values != None:
		for value in values.split(", "):
			if not os.path.isfile(value):
				raise optparse.OptionValueError("%s doesn't exist."%value)
	setattr(parser.values, option.dest, values)

def check_PERCENT (option, opt_str, value, parser):
	if not (0 <= value and value <= 100):
		 raise optparse.OptionValueError("%s option must be between 0 and 100."%opt_str)
	setattr(parser.values, option.dest, value)

def check_PROB (option, opt_str, value, parser):
	if not (0 <= value and value <= 1):
		 raise optparse.OptionValueError("%s option must be between 0 and 1."%opt_str)
	setattr(parser.values, option.dest, value)

def check_POS (option, opt_str, value, parser):
	if not (0 < value):
		 raise optparse.OptionValueError("%s option must be greater than 0."%opt_str)
	setattr(parser.values, option.dest, value)

def check_POSINTEGER (option, opt_str, value, parser):
	if not (0 < value and isinstance(value, (int, long))):
		 raise optparse.OptionValueError("%s option must be an integer greater than 0."%opt_str)
	setattr(parser.values, option.dest, value)

def check_NONNEG (option, opt_str, value, parser):
	if not (0 <= value):
		 raise optparse.OptionValueError("%s option must be greater than or equal to 0."%opt_str)
	setattr(parser.values, option.dest, value)

def check_NONNEGINTEGER (option, opt_str, value, parser):
	if not (0 <= value and isinstance(value, (int, long))):
		 raise optparse.OptionValueError("%s option must be an integer greater than or equal to 0."%opt_str)
	setattr(parser.values, option.dest, value)

def check_NONPOS (option, opt_str, value, parser):
	if not (value <= 0):
		 raise optparse.OptionValueError("%s option must be less than or equal to 0."%opt_str)
	setattr(parser.values, option.dest, value)

def check_N_AlignmentMerging (option, opt_str, value, parser):
	if not (-1 <= value and isinstance(value, (int, long))):
		 raise optparse.OptionValueError("%s option must be less than an integer greater than -1."%opt_str)
	setattr(parser.values, option.dest, value)

	
