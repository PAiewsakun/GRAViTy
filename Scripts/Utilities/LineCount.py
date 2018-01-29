def LineCount(filename):
	lines = 0
	buf_size = 1024 * 1024
	f = open(filename)
	read_f = f.read # loop optimization
	
	buf = read_f(buf_size)
	while buf:
		lines += buf.count('\n')
		buf = read_f(buf_size)
	
	return lines