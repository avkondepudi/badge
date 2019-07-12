import random
import argparse

def random_color():

	r = lambda: random.randint(0,255)
	return "#{:02x}{:02x}{:02x}".format(r(), r(), r())

def str2bool(x):
	if isinstance(x, bool): return x

	if x.lower() in ["true"]: return True
	elif x.lower() in ["false"]: return False
	else: raise argparse.ArgumentTypeError("Boolean value expected.")