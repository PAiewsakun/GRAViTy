def str2bool(v):
	if isinstance(v, bool):
		return v
	elif isinstance(v, str):
		return v.lower() in ("yes", "true", "t", "1")

