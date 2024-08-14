#! /bin/sh
import os, sys, re, string, glob
env = Environment(CCFLAGS = ['-O2'])  
libs = ['su','par','cwp','m'];
libpath = ['/usr/local/bin', '/bin','/usr/bin', '/usr/local/SU/lib']
cpppath = ['/usr/local/SU/include']

env.Program('acoustic_modeling_2d', 'acoustic_modeling_2d.c', LIBS=libs, LIBPATH=libpath,CPPPATH=cpppath)

