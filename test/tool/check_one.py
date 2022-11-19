#!/usr/bin/env python
import CompareFile as comf
import os
import sys
def red(txt):
    txt="\033[31m"+txt+"\033[0m"
    return txt
def green(txt):
    txt="\033[32m"+txt+"\033[0m"
    return txt
def yellow(txt):
    txt="\033[33m"+txt+"\033[0m"
    return txt
def check():
    testpass = comf.comparefile("result.txt", "result.ref", 1e-6)
    return testpass
def run():
    ierr=os.system("../../tool.exe < input > result.txt")
    if ierr == 0:
        return True
    else:
        return False

if __name__ == "__main__":
    # Check
    testpass = True
    testpass = testpass and run()
    testpass = testpass and check()
    if testpass:
        print("  One Processor      "+green( "PASS"))
    else:
        print("  One Processor      "+red("FAIL"))
