import signal

def timeout_handler(signum, frame):
	raise Exception("No input detected. Continue the program.")

def raw_input_with_timeout(prompt, timelimit = 5, default_input = "Y"):
	signal.signal(signal.SIGALRM, timeout_handler)
	signal.alarm(timelimit)
	input = default_input
	try:
		input = raw_input(prompt)
		signal.alarm(0)
	except Exception, exc: 
		print exc
	
	return input

