import time
import sys

def printRatio(i,n):
    print str(100*i/n)+'%',
    sys.stdout.flush()
    print "\r",
