#!/usr/bin/env pyhton
# encoding: utf-8

import os
path=input("please input a dir:\n")
print path
for root, dirs, files in os.walk(path):
    print "root:{},\ndirs:{},\nfiles:{}".format(root, dirs, files)
